import gzip
import heapq
import math
import statistics
import time
from collections import defaultdict, deque

ARQUIVO = "roadNet-PA.txt"  # ou 
N_EXEC = 30
Z_95 = 1.96  # Distribuição Normal (n >= 30, IC 95%)


def carregar_grafo(caminho):
    adj = defaultdict(list)
    arestas = []
    opener = gzip.open if caminho.endswith(".gz") else open
    with opener(caminho, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            u, v = map(int, line.split())
            adj[u].append(v)
            if u < v:
                arestas.append((u, v))
    return adj, arestas


def extrair_subgrafo(adj, inicio, limite):
    sub_nos = []
    vistos = {inicio}
    fila = deque([inicio])
    while fila and len(sub_nos) < limite:
        u = fila.popleft()
        sub_nos.append(u)
        for v in adj[u]:
            if v not in vistos:
                vistos.add(v)
                fila.append(v)
    mapa = {no: i for i, no in enumerate(sub_nos)}
    sub_adj = [[] for _ in range(len(sub_nos))]
    sub_arestas = []
    for u in sub_nos:
        iu = mapa[u]
        for v in adj[u]:
            if v in mapa:
                iv = mapa[v]
                sub_adj[iu].append(iv)
                if iu < iv:
                    sub_arestas.append((iu, iv))
    return sub_adj, sub_arestas


def benchmark(nome, func, n=N_EXEC):
    tempos = []
    resultado = None
    for _ in range(n):
        t0 = time.perf_counter()
        resultado = func()
        t1 = time.perf_counter()
        tempos.append(t1 - t0)
    media = statistics.mean(tempos)
    desvio = statistics.stdev(tempos)
    margem = Z_95 * (desvio / math.sqrt(n))
    print(f"{nome:<28} | Média: {media:.4f}s | DP: {desvio:.4f}s | IC 95%: [{media - margem:.4f}, {media + margem:.4f}]s")
    return resultado


# BFS
def bfs(adj, inicio):
    visitados = {inicio}
    fila = deque([inicio])
    while fila:
        u = fila.popleft()
        for v in adj[u]:
            if v not in visitados:
                visitados.add(v)
                fila.append(v)
    return len(visitados)


# DFS 
def dfs(adj, inicio):
    visitados = set()
    pilha = [inicio]
    while pilha:
        u = pilha.pop()
        if u not in visitados:
            visitados.add(u)
            for v in adj[u]:
                if v not in visitados:
                    pilha.append(v)
    return len(visitados)


# Verificação de Eulerianidade
def verificar_eulerianidade(adj):
    impares = sum(1 for viz in adj.values() if len(viz) % 2 != 0)
    inicio = next(iter(adj))
    alcancados = bfs(adj, inicio)
    conexo = (alcancados == len(adj))
    if not conexo:
        return f"Não euleriano (desconexo; {impares} vértices de grau ímpar)"
    if impares == 0:
        return "Circuito Euleriano"
    if impares == 2:
        return "Caminho Euleriano"
    return f"Não euleriano ({impares} vértices de grau ímpar)"


# Dijkstra (peso unitário = 1)
def dijkstra(adj, inicio):
    dist = {inicio: 0}
    heap = [(0, inicio)]
    while heap:
        d, u = heapq.heappop(heap)
        if d > dist[u]:
            continue
        nd = d + 1
        for v in adj[u]:
            if v not in dist or nd < dist[v]:
                dist[v] = nd
                heapq.heappush(heap, (nd, v))
    return len(dist)


# Bellman-Ford ( peso = 1)
def bellman_ford(n_nos, arestas, inicio=0):
    INF = 10**9
    dist = [INF] * n_nos
    dist[inicio] = 0
    for _ in range(n_nos - 1):
        atualizou = False
        for u, v in arestas:
            if dist[u] + 1 < dist[v]:
                dist[v] = dist[u] + 1
                atualizou = True
            if dist[v] + 1 < dist[u]:
                dist[u] = dist[v] + 1
                atualizou = True
        if not atualizou:
            break
    return max(d for d in dist if d < INF)


# Floyd-Warshall (no subgrafo reduzido, peso = 1)
def floyd_warshall(n_nos, arestas):
    INF = 10**9
    dist = [[INF] * n_nos for _ in range(n_nos)]
    for i in range(n_nos):
        dist[i][i] = 0
    for u, v in arestas:
        dist[u][v] = 1
        dist[v][u] = 1
    for k in range(n_nos):
        dk = dist[k]
        for i in range(n_nos):
            di = dist[i]
            dik = di[k]
            if dik == INF:
                continue
            for j in range(n_nos):
                if dik + dk[j] < di[j]:
                    di[j] = dik + dk[j]
    return dist[0][n_nos - 1]


# Algoritmo de Tarjan iterativo 
def tarjan_pontes(adj):
    disc = {}
    low = {}
    timer = 0
    pontes = 0

    for root in adj:
        if root in disc:
            continue
        disc[root] = low[root] = timer
        timer += 1
        pilha = [(root, -1, 0)]
        while pilha:
            u, pai, idx = pilha[-1]
            vizinhos = adj[u]
            if idx < len(vizinhos):
                pilha[-1] = (u, pai, idx + 1)
                v = vizinhos[idx]
                if v == pai:
                    continue
                if v not in disc:
                    disc[v] = low[v] = timer
                    timer += 1
                    pilha.append((v, u, 0))
                else:
                    if disc[v] < low[u]:
                        low[u] = disc[v]
            else:
                pilha.pop()
                if pai != -1:
                    if low[u] < low[pai]:
                        low[pai] = low[u]
                    if low[u] > disc[pai]:
                        pontes += 1
    return pontes


# Árvore Geradora Mínima - Kruskal com Union-Find (peso unitário = 1)
def kruskal(nos, arestas):
    pai = {u: u for u in nos}
    rank = {u: 0 for u in nos}

    def find(x):
        while pai[x] != x:
            pai[x] = pai[pai[x]]
            x = pai[x]
        return x

    arestas_agm = 0
    for u, v in arestas:
        ru, rv = find(u), find(v)
        if ru != rv:
            if rank[ru] < rank[rv]:
                pai[ru] = rv
            elif rank[ru] > rank[rv]:
                pai[rv] = ru
            else:
                pai[rv] = ru
                rank[ru] += 1
            arestas_agm += 1
    return arestas_agm


if __name__ == "__main__":
    print("Carregando grafo...")
    adj, arestas = carregar_grafo(ARQUIVO)
    nos = list(adj.keys())
    print(f"Carregado: {len(nos)} vértices, {len(arestas)} arestas.\n")

    sub500_adj, sub500_arestas = extrair_subgrafo(adj, 0, 500)
    sub200_adj, sub200_arestas = extrair_subgrafo(adj, 0, 200)

    print(f"Iniciando benchmarks (n = {N_EXEC} execuções, IC 95% via Normal Z = {Z_95}):\n")

    r_bfs = benchmark("BFS (grafo completo)", lambda: bfs(adj, 0))
    r_dfs = benchmark("DFS (grafo completo)", lambda: dfs(adj, 0))
    r_eul = benchmark("Eulerianidade (completo)", lambda: verificar_eulerianidade(adj))
    r_tar = benchmark("Tarjan Pontes (completo)", lambda: tarjan_pontes(adj))
    r_kru = benchmark("Kruskal AGM (completo)", lambda: kruskal(nos, arestas))
    r_dij = benchmark("Dijkstra (completo, w=1)", lambda: dijkstra(adj, 0))
    r_bf  = benchmark("Bellman-Ford (sub 500 nós)", lambda: bellman_ford(len(sub500_adj), sub500_arestas, 0))
    r_fw  = benchmark("Floyd-Warshall (sub 200 nós)", lambda: floyd_warshall(len(sub200_adj), sub200_arestas))

    print("\n--- Resultados obtidos ---")
    print(f"BFS alcançou:            {r_bfs} vértices (a partir do nó 0)")
    print(f"DFS alcançou:            {r_dfs} vértices (a partir do nó 0)")
    print(f"Eulerianidade:           {r_eul}")
    print(f"Tarjan (pontes):         {r_tar} pontes encontradas no grafo")
    print(f"Kruskal (floresta AGM):  {r_kru} arestas (|V| - 206 componentes)")
    print(f"Dijkstra alcançou:       {r_dij} vértices (a partir do nó 0)")
    print(f"Bellman-Ford (sub 500):  maior distância = {r_bf}")
    print(f"Floyd-Warshall (sub200): dist(0 -> 199)  = {r_fw}")