import gzip
import random
from collections import Counter, defaultdict, deque
import matplotlib.pyplot as plt
import networkx as nx

ARQUIVO = "roadNet-PA.txt"  # funciona também com "roadNet-PA.txt.gz"

# Leitura eficiente em lista de adjacência (grafo não-direcionado)
adj = defaultdict(list)
opener = gzip.open if ARQUIVO.endswith(".gz") else open

with opener(ARQUIVO, "rt", encoding="utf-8") as f:
    for line in f:
        if line.startswith("#"):
            continue
        u, v = map(int, line.split())
        adj[u].append(v)

n_vertices = len(adj)
graus = {u: len(viz) for u, viz in adj.items()}
n_arestas = sum(graus.values()) // 2

print(f"Vértices: {n_vertices}")
print(f"Arestas:  {n_arestas}")
print(f"Grau mínimo: {min(graus.values())}")
print(f"Grau máximo: {max(graus.values())}")
print(f"Grau médio:  {sum(graus.values()) / n_vertices:.4f}")
print(f"Densidade:   {(2 * n_arestas) / (n_vertices * (n_vertices - 1)):.6e}")

# Número de triângulos e Coeficiente de clusterização médio
adj_set = {u: set(viz) for u, viz in adj.items()}

soma_triangulos_por_no = 0
soma_clusterizacao = 0.0

for u, vizinhos in adj_set.items():
    k = len(vizinhos)
    if k < 2:
        continue
    viz_lista = list(vizinhos)
    tri_u = 0
    for i in range(k):
        v = viz_lista[i]
        viz_v = adj_set[v]
        for j in range(i + 1, k):
            if viz_lista[j] in viz_v:
                tri_u += 1
    soma_triangulos_por_no += tri_u
    soma_clusterizacao += (2 * tri_u) / (k * (k - 1))

# Cada triângulo é contado 3 vezes (uma em cada vértice)
n_triangulos = soma_triangulos_por_no // 3
clusterizacao_media = soma_clusterizacao / n_vertices

print(f"\nNúmero de triângulos:               {n_triangulos}")
print(f"Coeficiente de clusterização médio: {clusterizacao_media:.4f}")

# Distribuição de graus
dist_graus = Counter(graus.values())
print("\nDistribuição de graus (grau: quantidade):")
for g in sorted(dist_graus):
    print(f"  Grau {g}: {dist_graus[g]} ({dist_graus[g]/n_vertices*100:.2f}%)")


# Componentes conexas 
visitados = set()
componentes = []

for node in adj:
    if node not in visitados:
        comp = []
        fila = deque([node])
        visitados.add(node)
        while fila:
            atual = fila.popleft()
            comp.append(atual)
            for viz in adj[atual]:
                if viz not in visitados:
                    visitados.add(viz)
                    fila.append(viz)
        componentes.append(comp)

componentes.sort(key=len, reverse=True)
print(f"\nNúmero de componentes conexas: {len(componentes)}")
print("Tamanho das 10 maiores componentes:", [len(c) for c in componentes[:10]])
print("Resumo de tamanhos (tamanho: qtd de componentes):", Counter(len(c) for c in componentes))

# Estimativa de Raio e Comprimento Médio dos Caminhos 
maior_comp = componentes[0]
amostra = random.sample(maior_comp, 100) 

excentricidades = []
soma_dist = 0
total_pares = 0

for origem in amostra:
    dist = {origem: 0}
    fila = deque([origem])
    max_d = 0
    while fila:
        atual = fila.popleft()
        d_atual = dist[atual]
        if d_atual > max_d:
            max_d = d_atual
        for viz in adj[atual]:
            if viz not in dist:
                dist[viz] = d_atual + 1
                fila.append(viz)
    excentricidades.append(max_d)
    soma_dist += sum(dist.values())
    total_pares += len(dist) - 1

print(f"\nRaio estimado (mín. excentricidade na amostra): {min(excentricidades)}")
print(f"Comprimento médio dos caminhos estimado:        {soma_dist / total_pares:.2f}")

# Gráfico da distribuição de graus
plt.figure(figsize=(8, 4))
xs = sorted(dist_graus.keys())
ys = [dist_graus[x] for x in xs]
plt.bar(xs, ys, color="#2b6cb0")
plt.xlabel("Grau")
plt.ylabel("Número de vértices")
plt.title("Distribuição de Graus - roadNet-PA")
plt.xticks(xs)
plt.tight_layout()
plt.savefig("distribuicao_graus.png", dpi=150)

# Visualização reduzida (subgrafo BFS de 300 nós a partir do nó 0)
sub_nos = []
fila = deque([0])
vistos = {0}
while fila and len(sub_nos) < 300:
    u = fila.popleft()
    sub_nos.append(u)
    for v in adj[u]:
        if v not in vistos:
            vistos.add(v)
            fila.append(v)

sub_set = set(sub_nos)
G_sub = nx.Graph()
for u in sub_nos:
    for v in adj[u]:
        if v in sub_set:
            G_sub.add_edge(u, v)

plt.figure(figsize=(8, 8))
pos = nx.spring_layout(G_sub, seed=42)
nx.draw(G_sub, pos, node_size=15, node_color="#e53e3e", edge_color="#a0aec0", width=0.6)
plt.title("Visualização Reduzida (300 vértices via BFS) - roadNet-PA")
plt.tight_layout()
plt.savefig("visualizacao_roadNet_PA.png", dpi=150)
print("\nGráficos salvos: distribuicao_graus.png e visualizacao_roadNet_PA.png")