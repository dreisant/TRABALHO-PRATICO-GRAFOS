# Análise de Grafos do SNAP — Pennsylvania Road Network (`roadNet-PA`)

Trabalho Prático da disciplina de **Teoria dos Grafos** — Instituto de Computação (IC), Universidade Federal da Bahia (UFBA).

Este repositório contém os scripts em Python, o relatório em LaTeX (formato WebMedia) e os gráficos gerados na análise estrutural, algorítmica e de robustez da malha rodoviária do estado da Pensilvânia (`roadNet-PA`), disponibilizada pelo *Stanford Network Analysis Project* (SNAP).

---

## Sobre o Conjunto de Dados

* **Fonte:** [SNAP — Road network of Pennsylvania (`roadNet-PA`)](https://snap.stanford.edu/data/roadNet-PA.html)
* **Modelagem:** Vértices representam cruzamentos e extremidades de vias; arestas não-direcionadas representam os trechos de estradas que os conectam.
* **Dimensões:** `1.088.092` vértices e `1.541.898` arestas não-direcionadas.

### Tratamento dos Dados
* **Grafo não-direcionado simples:** o arquivo bruto registra `3.083.796` linhas por listar cada aresta nos dois sentidos ($u \to v$ e $v \to u$). Os scripts consolidam a leitura em lista de adjacência como grafo não-direcionado simples.
* **Lista de adjacência:** adotada em memória devido à extrema esparsidade ($D \approx 2,60 \times 10^{-6}$).
* **Maior componente conexa (GCC):** extraída para estimativas de distância e simulações de robustez (`1.087.562` vértices — 99,95% da rede).
* **Subgrafos reduzidos (BFS):** extraídos a partir do vértice `0` para visualização (`300` nós) e para algoritmos computacionalmente inviáveis no grafo completo (Bellman-Ford em `500` nós; Floyd-Warshall em `200` nós).
* **Pesos unitários ($w = 1$):** atribuídos às arestas nos algoritmos ponderados (Dijkstra, Bellman-Ford, Floyd-Warshall e Kruskal).

---

## Estrutura do Repositório

```text
.
├── analise.py                   # Parte I: análise estrutural, graus, componentes, triângulos e visualização
├── parte2.py                    # Parte II: 8 algoritmos clássicos + benchmarks (n=30, IC 95%)
├── parte3.py                    # Parte III: Small-World, Lei de Potência (log-log) e Robustez (5%)
├── roadNet-PA.txt.gz            # Conjunto de dados do SNAP (ou roadNet-PA.txt)
├── distribuicao_graus.png       # Gráfico gerado pela Parte I
├── visualizacao_roadNet_PA.png  # Subgrafo reduzido (300 nós) gerado pela Parte I
├── power_law_loglog.png         # Gráfico log-log gerado pela Parte III
├── robustez_comparacao.png      # Boxplots de robustez gerados pela Parte III
└── README.md                    # Este arquivo
```

---

## Pré-requisitos e Instalação

1. **Python 3.8+** instalado no sistema.
2. Clone este repositório (ou baixe os arquivos) e acesse a pasta pelo terminal:
   ```bash
   git clone https://github.com/SEU-USUARIO/SEU-REPOSITORIO.git
   cd SEU-REPOSITORIO
   ```
3. Instale as bibliotecas externas utilizadas para geração de gráficos e visualização:
   ```bash
   pip install matplotlib networkx
   ```
   *(Observação: todos os algoritmos da Parte II e os cálculos estruturais das Partes I e III foram implementados com estruturas nativas do Python, sem depender de bibliotecas prontas de grafos. O `networkx` é usado apenas para o layout da visualização reduzida na Parte I).*

4. **Arquivo de dados:** certifique-se de que `roadNet-PA.txt` (ou `roadNet-PA.txt.gz`) esteja na raiz do projeto. Caso precise baixá-lo novamente, acesse [snap.stanford.edu/data/roadNet-PA.html](https://snap.stanford.edu/data/roadNet-PA.html).
   * Se utilizar o arquivo compactado `.txt.gz`, ajuste a constante `ARQUIVO = "roadNet-PA.txt.gz"` no início de cada script.

---

## Como Executar

Execute cada script separadamente na raiz do projeto:

### Parte I — Análise Estrutural Obrigatória
Calcula vértices, arestas, graus, distribuição de graus, densidade, componentes conexas, triângulos, coeficiente de clusterização médio, estimativas de raio e caminho médio (via BFS) e gera os gráficos iniciais.

```bash
python analise.py
```
* **Tempo estimado:** ~1 a 2 minutos.
* **Saídas geradas:** estatísticas no terminal + `distribuicao_graus.png` e `visualizacao_roadNet_PA.png`.

### Parte II — Algoritmos da Disciplina e Benchmark
Executa `n = 30` repetições de BFS, DFS iterativa, Verificação de Eulerianidade, Tarjan (pontes), Kruskal (AGM), Dijkstra, Bellman-Ford (subgrafo de 500 nós) e Floyd-Warshall (subgrafo de 200 nós), reportando média, desvio padrão e Intervalo de Confiança de 95% (distribuição Normal, $Z = 1,96$).

```bash
python parte2.py
```
* **Tempo estimado:** ~3 a 5 minutos.
* **Saídas geradas:** tabela de tempos e resultados algorítmicos no terminal.

### Parte III — Small-World, Lei de Potência e Robustez
Avalia a propriedade *small-world* (sem biblioteca pronta), plota a distribuição de graus em escala log-log e simula a remoção de 5% dos vértices da GCC (`54.378` nós) em $T = 30$ repetições aleatórias comparadas ao ataque direcionado por centralidade de grau.

```bash
python parte3.py
```
* **Tempo estimado:** ~1 a 2 minutos.
* **Saídas geradas:** métricas no terminal + `power_law_loglog.png` e `robustez_comparacao.png`.

---

## Resumo dos Principais Resultados

| Propriedade / Análise | Resultado |
| :--- | :--- |
| **Vértices / Arestas** | `1.088.092` / `1.541.898` |
| **Grau (mín / máx / médio)** | `1` / `9` / `2,8341` |
| **Densidade** | `2,6047 × 10⁻⁶` |
| **Componentes conexas** | `206` (GCC com `1.087.562` nós — 99,95%) |
| **Diâmetro / Raio / Caminho médio** | `786` / `≤ 417` / `≈ 309,35` |
| **Clusterização média / Triângulos** | `0,0465` / `67.150` |
| **Eulerianidade** | Não euleriano (`728.846` vértices de grau ímpar) |
| **Pontes (Tarjan)** | `216.994` arestas de corte (`14,07%`) |
| **Small-World / Lei de Potência** | Não / Não |
| **Robustez (5% aleatório vs. direcionado)** | GCC retém `98,18%` (aleatório) vs. `86,26%` (direcionado) |

---

## Autor

**Davi Reis Santos Silva**  
Universidade Federal da Bahia (UFBA) — Instituto de Computação