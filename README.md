# ABC_2026

[Lucas Faria](https://github.com/lcsfaria), [Camila Sinimbú Forte](https://github.com/camiliane), [Carolina Silva de Carvalho](https://github.com/carolinascarvalho), [Glen Garcia](https://github.com/glenjasper), [Luca Vasconcelos da Gama Garcia](https://github.com/lucagamam), [Ivan Gomes Cruz](https://github.com/IvanGomesCruz).

EM CONSTRUÇÃO 🏗️🏗️🏗️🏗️🏗️🏗️🏗️🏗️

## Introdução 

### Computação Bayesiana Aproximada (ABC)
O ABC foi introduzido por Tavaré et al. (1997) e Pritchard et al. (1999) para resolver problemas de likelihood intratáveis na genética de populações. O método aborda o problema de encontrar uma distribuição posterior dos parâmetros em um modelo que explica um conjunto de dados complexo. 

O princípio fundamental é: Dada uma amostra de parâmetros $θ$, um conjunto de dados $D'$ é simulado em função deles. Se o conjunto de  $D'$ é muito diferente do dado observado $D$, então os $θ$ responsáveis por $D'$ são descartados. Dada uma função de distância $p$ e uma tolerância ε, $D'$ é aceito caso $ρ(D,D') ≤ ε$. Por fim, após a seleção dos cenários que tiveram o menor $ρ(D,D')$, volta-se para os valores dos nossos parâmetros $θ$ que constituem nossas distribuições a posteriori (Figura 1).

<img width="1744" height="675" alt="image" src="https://github.com/user-attachments/assets/9fc1486c-b952-4677-905a-99bd34b74f9f" />

Figura 1. **Estimação de parâmetros através do ABC**. Fluxo da estimativa de parâmetros a partir de dados observados e simulados no ABC. 

### Inferência da Dinâmica da Miscigenação 
Nossa análise é estruturada em quatro etapas: i) geração dos cenários demográficos candidatos; ii) simulação do modelo para cada cenário candidato; iii) construção das estatísticas sumárias dos dados simulados; e iv) construção da posteriori, dado o cálculo das distâncias entre as estatísticas sumárias simuladas e observadas e a retenção dos cenários mais próximos.

<img width="617" height="675" alt="ABC" src="https://github.com/user-attachments/assets/6fe3cd67-76b1-4d10-8803-b94db53bbb86" />

---

## Fluxo e scripts:
O pipeline segue da seguinte forma:

**Geração de cenários → Simulação dos tracts → Construção dos histogramas → Comparação das distâncias**

- [generate_scenarios.py](#Geração-de-cenários): Geração dos cenários simulados
- [simulate_tracts.py](#Simulação-de-tracts): Simulação dos tracts de ancestralidade para cada cenário e réplica (Liang&Nielsen)
- [build_tract_bins.py](#Construção-dos-histogramas): Construção de histogramas por ancestralidade (somente para modos com bins)
- [compute_tract_distances.py](#Cálculo-das-distâncias-entre-simulado-e-observado): Cálculo das distâncias entre simulado e observado
- [plots.py](#Plots): Plots (gráficos de distribuições posteriores) e diagnósticos (KDE e HPD)

### Modos de execução
(Parâmetro: mode)
- **expanded**: Expande os histogramas em distribuições contínuas repetindo os bins conforme seus counts. Preserva a forma e a frequência dos tracts.
- **count_based**: Compara apenas os valores de counts dos bins, ignorando os comprimentos reais dos tracts. Serve principalmente para testes/comparação metodológica;
- **raw_tracts**: Usa diretamente os comprimentos brutos dos tracts simulados e observados, sem histogramas nem discretização. É o modo mais fiel aos dados originais.

### Modos de agregação entre as ancestralidades
(Parâmetro: aggregation)

Após calcular as 5 estatísticas (parâmetro: statistics) por ancestralidade, os valores são consolidados em 1 valor por cenário.

- **arithmetic** (média simples): Todas as ancestralidades têm o mesmo peso;
- **weighted** (ponderada por counts): Mais tracts = mais peso na métrica;
- **euclidean** (quadrática): Penaliza fortemente discrepâncias grandes (em qualquer ancestralidade)

### Arquivo config.ini
```
[Demography]
n_generations = 21
founder_population = NAT, AFR, EUR
generation_times = 18, 17, 16, 12, 11, 10, 6, 5, 4

[Populations]
ancestry_proportions = EUR 0.783, AFR 0.144, NAT 0.073

[Genome]
chromosomes = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22
n_chromosomes = 7472
n_bins = 50
window_size = 0.1
recombination_map = /resources/recombination.txt

[Simulation]
n_scenarios = 3000000
n_replicas = 1
n_threads = 40

[ObservedData]
reference_tracts = /LAI/observado #caminho para os arquivos de tracts observados
scenario_percentage = 1

[Metrics]
mode = expanded #modos de execução
statistics = wasserstein, ks_statistic, energy, jensen_shannon, hellinger #métricas de distância
aggregation = euclidean #agregação entre as ancestralidades
```
---
## Orquestrador do pipeline
**Para rodar todo o pipeline:**
```
python run_pipeline.py -c config.ini --seed 123 -o outputs
```

## Geração de cenários
**1º passo do pipeline**
```
python generate_scenarios.py -c config.ini --seed 123 -o outputs
```
Gera cenários candidatos utilizando uma amostragem de Dirichlet. Cada cenário representa um conjunto de pulsos migratórios provenientes de populações parentais, que ocorrem em gerações específicas.

Os cenários são filtrados com base nas proporções finais de ancestralidade que produzem, garantindo a compatibilidade com a composição de ancestralidade esperada definida no arquivo de configuração.

#### Inputs:
Do arquivo [config.ini](#Arquivo-config.ini)
```
[Demography]
n_generations = 21
founder_population = NAT, AFR, EUR
generation_times = 18, 17, 16, 12, 11, 10, 6, 5, 4

[Populations]
ancestry_proportions = EUR 0.783, AFR 0.144, NAT 0.073

[Genome]
window_size = 0.1

[Simulation]
n_scenarios = 3000000
n_threads = 40
```

#### Outputs:
- *informative_prior.txt*: Lista de comandos para a simulação dos tracts (Liang-Nielsen) correspondentes aos cenários aceitos.
- *all_scenarios.txt*: Lista bruta de todos os cenários amostrados, antes da filtragem.
- *informative_prior_onlyValues.txt*: Matriz contendo apenas as proporções de migração (valores m) dos cenários aceitos, ordenadas por população e pulso (da mais recente para a mais antiga). Este arquivo é utilizado posteriormente durante a seleção a posteriori.
- *informative_prior_by_generations.txt*

<img width="1917" height="996" alt="image" src="https://github.com/user-attachments/assets/bd902003-67db-433d-aeca-b7a21e7b57bf" />

--- 

## Simulação de tracts
**2º passo do pipeline**
```
python simulate_tracts.py -c config.ini -s outputs/informative_prior.txt --seed 123 -o outputs/simulated_tracts
```
Simulação dos tracts de ancestralidade ao longo dos cromossomos para cada cenário demográfico aceito usando um simulador de miscigenação baseado em Wright-Fisher (Liang-Nielsen). 

Cada simulação produz uma lista de tracts de ancestralidade com seus respectivos comprimentos.

#### Inputs:
- Do arquivo [config.ini](#Arquivo-config.ini)
```
[Genome]
chromosomes = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22
recombination_map = /resources/recombination.txt

[Simulation]
n_scenarios = 3000000
n_replicas = 1
n_threads = 40
```
- Arquivo de cenários: *informative_prior.txt*.

#### Outputs:
Estrutura do diretório:
```
outputs/
  simulated_tracts/
    Scenario_1/
        CHR1.rep1.tracts.txt
        CHR1.rep2.tracts.txt
        ...
        CHR2.rep1.tracts.txt
        ...
    Scenario_2/
        ...
```
Cada arquivo contém os tracts de ancestralidade simulada com o seguinte formato:

    population   tract_length   chromosome_copy   chromosome_id

---

## Construção dos histogramas
> **Obs**: Nos modos expanded ou count_based


**3º passo do pipeline**
```
python build_tract_bins.py -c config.ini -b outputs/simulated_tracts -o outputs/tract_bins
```
Conversão dos tracts de ancestralidade simulados em intervalos de histograma de comprimento fixo. Para cada cromossomo simulado, os comprimentos dos trechos são agrupados em intervalos de acordo com o número de intervalos especificado no arquivo de configuração. 

Produz uma distribuição de comprimento de trecho para cada população ancestral. Esses histogramas são posteriormente usados ​​para comparar os dados simulados com as distribuições de trechos de referência observadas.

#### Inputs:
- Do arquivo [config.ini](#Arquivo-config.ini)
```
[Genome]
chromosomes = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22
n_bins = 50
recombination_map = /resources/recombination.txt

[Simulation]
n_threads = 40

[Metrics]
mode = expanded, count_based
```
- Caminho para os arquivos dos tracts simulados: outputs/simulated_tracts.
```
outputs/
  simulated_tracts/
      Scenario_X/
        CHR*.rep*.tracts.txt
```

#### Outputs:
Estrutura do diretório:
```
outputs/
  tract_bins/
    Scenario_X/
        CHR1.rep1.bins.txt
        CHR1.rep2.bins.txt
        ...
```
Cada arquivo contém:

    population    bin_length    tract_count
    
- Em que:
  - bin_length = limite superior do bin
  - tract_count = número de tracts que estão contidos nesse intervalo


<img width="1435" height="972" alt="image" src="https://github.com/user-attachments/assets/07095e4a-66eb-48bf-a7c7-262705ec3170" />

---

## Cálculo das distâncias entre simulado e observado
**Último passo do pipeline**

> **Obs**: Nos modos expanded ou count_based (baseado em bins)
```
python compute_tract_distances.py -c config.ini --tracts_dir outputs/tract_bins --prior_file outputs/informative_prior_onlyValues.txt --prior_generations_file informative_prior_by_generations.txt -o outputs/tract_distances
```

**ou**


> **Obs**: No modo raw_tract (baseado em tracts)
```
python compute_tract_distances.py -c config.ini --tracts_dir outputs/simulated_tracts --prior_file outputs/informative_prior_onlyValues.txt --prior_generations_file informative_prior_by_generations.txt -o outputs/tract_distances
```

Comparação das distribuições dos tracts simulados com as distribuições de referência observadas usando múltiplas métricas de distância estatística. Para cada réplica da simulação e cromossomo, as distribuições de comprimento das sequências são comparadas aos dados observados. Os cenários são classificados de acordo com sua distância em relação à distribuição observada. 

Esta etapa implementa o estágio de seleção por Computação Bayesiana Aproximada (ABC) do pipeline.

#### Inputs:
- Do arquivo [config.ini](#Arquivo-config.ini)
```
[Simulation]
n_threads = 40

[ObservedData]
reference_tracts = /LAI/observado #caminho para os arquivos dos tracts observados 
scenario_percentage = 1

[Metrics]
statistics = wasserstein, ks_statistic, energy, jensen_shannon, hellinger #métricas de distância
aggregation = euclidean #agregação entre as ancestralidades
```

- Caminho para os arquivos dos tracts simulados ou dos bins: outputs/simulated_tracts **ou** outputs/tract_bins.
```
outputs/
  tract_bins/
    Scenario_X/
       CHR*.rep*.<bins|tracts>.txt
```

- Arquivo de cenários contendo apenas as proporções de migração (valores m): *informative_prior_onlyValues.txt*.
- Arquivo contendo os cenários por geração: *informative_prior_by_generations.txt*

#### Outputs:
- *all_metrics.txt*: Tabela contendo todas as métricas de distância calculadas para cada cenário, cromossomo e réplica.
- **ranking_CHR*_rep*_metric.csv: Arquivos de classificação (por cromossomo e réplica).
- **posterior_CHR*_rep*_metric.txt: Arquivo das posteriores contendo os cenários com melhor classificação de acordo com a porcentagem de cenário especificada no arquivo de configuração.
- **posterior_generations_CHR*_rep*_<metric>.txt: Arquivo das posteriores contendo os cenários com melhor classificação por geração.


- Métricas calculadas: 
  - Distância de Wasserstein
  - Estatística de Kolmogorov-Smirnov
  - Distância de energia
  - Divergência de Jensen-Shannon
  - Diferença média
  - Diferença mediana
  - Diferença de desvio padrão
    
---

## Plots 
```
python plots.py -c config.ini --posterior_dir outputs/tract_distances --prior_file outputs/informative_prior_onlyValues.txt -o outputs/plots
```
Gera os gráficos de comparação entre a distribuição posterior e a distribuição a priori para os resultados do ABC. 

As distribuições posteriores são visualizadas usando a estimativa de densidade de kernel (KDE), e a incerteza é resumida usando intervalos de densidade posterior mais alta (HPD).

#### Inputs:
- Do arquivo [config.ini](#Arquivo-config.ini): Utilizado para recuperar nomes e estrutura da população.

- Caminho para os arquivos das posteriores
```
outputs/
  tract_distances/
      CHR*/
        posterior_CHR*_rep*_<metric>.txt
        posterior_generations_CHR*_rep*_<metric>.txt
```

- Arquivo de cenários contendo apenas as proporções de migração (valores m): *informative_prior_onlyValues.txt*.

#### Outputs:
Plots salvos como:
```
outputs/
    plots/
        posterior_*_kde.png
        posterior_*_hpd.png
        posterior_generations_*_hpd.png
```
  - Gráficos KDE (distribuições a posteriori versus a priori com intervalos HPD)
  - Gráficos de intervalo HPD ao longo dos pulsos de miscigenação
  - Gráficos de intervalo HPD ao longo das gerações de miscigenação

<img width="1733" height="743" alt="image" src="https://github.com/user-attachments/assets/576e54f9-089b-4f28-86cb-77e40fa7743b" />


<img width="1376" height="841" alt="image" src="https://github.com/user-attachments/assets/bbbf1f77-03d8-48d1-a9d4-0265a8476233" />

####Validação
EM CONSTRUÇÃO 🏗️🏗️🏗️🏗️🏗️🏗️🏗️🏗️
