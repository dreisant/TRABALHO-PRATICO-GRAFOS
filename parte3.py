import gzip
import math
import random
import statistics
from collections import Counter, defaultdict, deque
import matplotlib.pyplot as plt

ARQUIVO = "roadNet-PA.txt" 
SEED = 42
T_REP = 30  # repetições para remoção aleatória (Apêndice C)


def carregar_grafo(caminho):
    adj = defaultdict(list)
    opener = gzip.open if caminho.endswith(".gz") else open
    with opener(caminho, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            u, v = map(int, line.split())
            adj[u].append(v)
    return adj


def extrair_maior_componente(adj):
    visitados = set()
    maior = []
    for node in adj:
        if node not in visitados:
            comp = []
            fila = deque([node])
            visitados.add(node)
            while fila:
                u = fila.popleft()
                comp.append(u)
                for v in adj[u]:
                    if v not in visitados:
                        visitados.add(v)
                        fila.append(v)
            if len(comp) > len(maior):
                maior = comp
    maior_set = set(maior)
    adj_gcc = {u: [v for v in adj[u] if v in maior_set] for u in maior}
    return adj_gcc


def avaliar_rede(adj_gcc, removidos, n_bfs_dist=5):
    """Calcula Métricas A, B, C e D após remoção de um conjunto de vértices."""
    visitados = set(removidos)
    n_restantes = len(adj_gcc) - len(removidos)
    maior_comp = []
    n_componentes = 0
    isolados = 0

    for u in adj_gcc:
        if u in visitados:
            continue
        tem_vizinho = False
        for v in adj_gcc[u]:
            if v not in removidos:
                tem_vizinho = True
                break
        if not tem_vizinho:
            visitados.add(u)
            n_componentes += 1
            isolados += 1
            if not maior_comp:
                maior_comp = [u]
            continue
        comp = []
        fila = deque([u])
        visitados.add(u)
        while fila:
            atual = fila.popleft()
            comp.append(atual)
            for v in adj_gcc[atual]:
                if v not in visitados:
                    visitados.add(v)
                    fila.append(v)
        n_componentes += 1
        if len(comp) > len(maior_comp):
            maior_comp = comp

    amostra = random.sample(maior_comp, min(n_bfs_dist, len(maior_comp)))
    soma_dist = 0
    total_pares = 0
    for origem in amostra:
        dist = {origem: 0}
        fila = deque([origem])
        while fila:
            atual = fila.popleft()
            d_atual = dist[atual]
            for v in adj_gcc[atual]:
                if v not in removidos and v not in dist:
                    dist[v] = d_atual + 1
                    fila.append(v)
        soma_dist += sum(dist.values())
        total_pares += len(dist) - 1

    dist_media = soma_dist / total_pares if total_pares > 0 else 0.0
    frac_maior = len(maior_comp) / n_restantes
    frac_isolados = isolados / n_restantes
    return len(maior_comp), frac_maior, n_componentes, dist_media, frac_isolados


if __name__ == "__main__":
    random.seed(SEED)
    print("Carregando grafo...")
    adj = carregar_grafo(ARQUIVO)
    n_total = len(adj)
    graus = {u: len(viz) for u, viz in adj.items()}
    k_medio = sum(graus.values()) / n_total

    # ITEM 1: SMALL-WORLD 
    L_obs = 309.35   # obtido na Parte I (500 amostras BFS)
    C_obs = 0.0465   # obtido na Parte I (cálculo exato)
    L_rand = math.log(n_total) / math.log(k_medio)
    C_rand = k_medio / n_total
    sigma = (C_obs / C_rand) / (L_obs / L_rand)

    print("\n=== 1. SMALL-WORLD ===")
    print(f"L(G) observado:  {L_obs:.2f}  | L_random (ln N / ln k): {L_rand:.2f}  | Razão L/L_rand: {L_obs/L_rand:.2f}x")
    print(f"C(G) observado:  {C_obs:.4f} | C_random (k / N):       {C_rand:.6e} | Razão C/C_rand: {C_obs/C_rand:.1f}x")
    print(f"Coeficiente sigma = (C/C_rand)/(L/L_rand): {sigma:.2f}")

    # ITEM 2: LEI DE POTÊNCIA 
    dist_graus = Counter(graus.values())
    ks = sorted(dist_graus.keys())
    ps = [dist_graus[k] / n_total for k in ks]

    plt.figure(figsize=(7, 4.5))
    plt.loglog(ks, ps, "o-", color="#2b6cb0", markersize=6, linewidth=1.5)
    plt.xlabel("Grau k (escala log)")
    plt.ylabel("Probabilidade P(k) (escala log)")
    plt.title("Distribuição de Graus em Escala Log-Log - roadNet-PA")
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig("power_law_loglog.png", dpi=150)
    print("\n=== 2. LEI DE POTÊNCIA ===")
    print("Gráfico salvo: power_law_loglog.png")

    # ITEM 3: ROBUSTEZ 

    print("\n=== 3. ROBUSTEZ (extraindo maior componente conexa...) ===")
    adj_gcc = extrair_maior_componente(adj)
    nos_gcc = list(adj_gcc.keys())
    n_gcc = len(nos_gcc)
    r = int(0.05 * n_gcc)
    print(f"GCC: {n_gcc} vértices | Removendo 5% = {r} vértices")

    # 3a. Remoção aleatória (T = 30 repetições)
    print(f"Executando {T_REP} repetições de remoção aleatória de 5%...")
    rand_A_abs, rand_A_frac, rand_B, rand_C, rand_D = [], [], [], [], []

    for rep in range(T_REP):
        removidos = set(random.sample(nos_gcc, r))
        tam_m, frac_m, n_comp, d_med, f_iso = avaliar_rede(adj_gcc, removidos, n_bfs_dist=5)
        rand_A_abs.append(tam_m)
        rand_A_frac.append(frac_m)
        rand_B.append(n_comp)
        rand_C.append(d_med)
        rand_D.append(f_iso)
        if (rep + 1) % 10 == 0:
            print(f"  Repetição {rep + 1}/{T_REP} concluída...")

    # 3b. remoção dos centrais
    print("Executando ataque direcionado (top 5% por centralidade de Grau)...")
    nos_ordenados = sorted(nos_gcc, key=lambda u: len(adj_gcc[u]), reverse=True)
    removidos_cent = set(nos_ordenados[:r])
    cent_A_abs, cent_A_frac, cent_B, cent_C, cent_D = avaliar_rede(adj_gcc, removidos_cent, n_bfs_dist=30)

    def resumo(nome, vals, val_cent, fmt="{:.4f}"):
        m = statistics.mean(vals)
        dp = statistics.stdev(vals)
        print(f"  {nome:<28} | Aleatório: {fmt.format(m)} ± {fmt.format(dp)} | Direcionado: {fmt.format(val_cent)}")

    print("\n--- Comparação de Robustez (5% removidos na GCC) ---")
    resumo("Métrica A (Tam. maior comp.)", rand_A_abs, cent_A_abs, "{:.0f}")
    resumo("Métrica A (Fração S)",         rand_A_frac, cent_A_frac, "{:.4f}")
    resumo("Métrica B (Nº componentes)",   rand_B, cent_B, "{:.1f}")
    resumo("Métrica C (Dist. média)",      rand_C, cent_C, "{:.2f}")
    resumo("Métrica D (Fração isolados)",  rand_D, cent_D, "{:.4f}")


    # Boxplots comparativos (Aleatório vs Direcionado)
    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    dados = [
        (axes[0, 0], rand_A_frac, cent_A_frac, "Métrica A: Fração da Maior Componente (S)"),
        (axes[0, 1], rand_B,      cent_B,      "Métrica B: Número de Componentes Conexas"),
        (axes[1, 0], rand_C,      cent_C,      "Métrica C: Comprimento Médio dos Caminhos"),
        (axes[1, 1], rand_D,      cent_D,      "Métrica D: Fração de Nós Isolados"),
    ]
    for ax, vals, v_cent, titulo in dados:
        ax.boxplot(vals, widths=0.4, patch_artist=True, boxprops=dict(facecolor="#90cdf4"))
        ax.axhline(v_cent, color="#e53e3e", linestyle="--", linewidth=2, label=f"Ataque direcionado ({v_cent:.4g})")
        ax.set_title(titulo, fontsize=10)
        ax.set_xticks([1])
        ax.set_xticklabels(["Remoção Aleatória (5%)"])
        ax.legend(fontsize=8)
        ax.grid(axis="y", linestyle="--", alpha=0.5)

    plt.suptitle("Análise de Robustez - roadNet-PA (Aleatória vs. Ataque por Grau)", fontsize=12)
    plt.tight_layout()
    plt.savefig("robustez_comparacao.png", dpi=150)
    print("\nGráfico salvo: robustez_comparacao.png")