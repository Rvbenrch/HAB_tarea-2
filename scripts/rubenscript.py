#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Propagación en redes con DIAMOnD (versión didáctica, sin dependencias raras).
- Carga red (STRING HUGO / DIAMOnD CSV / GUILD con pesos)
- Lee semillas (ENO1, PGK1, HK2 por defecto)
- Ejecuta DIAMOnD (hipergeométrica sobre vecinos del módulo)
- Exporta ranking, metadatos, figura y subred sugerida

Referencias:
- Ghiassian SD, Menche J, Barabási A-L. A DIseAse MOdule Detection (DIAMOnD) algorithm (2015).
Implementación simplificada para docencia.

Uso:
python scripts/tu_script.py --network data/string_network_filtered_hugo-400.tsv --network-format string \
    --seeds data/genes_seed.txt --k 200 --outdir results
"""

import argparse
import json
import math
import os
from datetime import datetime
from typing import Iterable, Tuple, Dict, Set, List

import pandas as pd
import networkx as nx
import numpy as np

# SciPy es opcional; si está, usamos su hipergeométrica
try:
    from scipy.stats import hypergeom
    _HAVE_SCIPY = True
except Exception:
    _HAVE_SCIPY = False


def comb(n: int, k: int) -> float:
    """Combinatoria segura (para fallback sin SciPy)."""
    if k < 0 or k > n:
        return 0.0
    return math.comb(n, k)


def hypergeom_sf(k: int, N: int, K: int, n: int) -> float:
    """
    Cola superior (P[X >= k]) de la hipergeométrica para DIAMOnD.
    N: tamaño población (nodos-1)
    K: 'éxitos' en la población (tamaño módulo)
    n: draws (grado del nodo)
    k: enlaces observados del nodo al módulo
    """
    if _HAVE_SCIPY:
        # survival function: P(X >= k)
        return float(hypergeom.sf(k - 1, N, K, n))
    # Fallback exacto por suma de combinatorias
    # P = sum_{i=k..min(n,K)} [C(K,i)*C(N-K, n-i)] / C(N,n)
    denom = comb(N, n)
    max_i = min(n, K)
    s = 0.0
    for i in range(k, max_i + 1):
        s += comb(K, i) * comb(N - K, n - i)
    return s / denom if denom > 0 else 1.0


# ------------------ Carga de redes ------------------

def load_string_network(path: str) -> nx.Graph:
    df = pd.read_csv(path, sep="\t")
    # Esperamos columnas protein1_hugo, protein2_hugo, combined_score
    ucol = [c for c in df.columns if "protein1" in c][0]
    vcol = [c for c in df.columns if "protein2" in c][0]
    G = nx.Graph()
    for _, r in df.iterrows():
        u, v = str(r[ucol]), str(r[vcol])
        if u == v:
            continue
        w = float(r.get("combined_score", 1.0))
        # normalizamos a [0,1] si viene como 0..1000
        if w > 1.0:
            w = w / 1000.0
        G.add_edge(u, v, weight=w)
    return G


def load_diamond_network(path: str) -> nx.Graph:
    # CSV: u,v (IDs)
    df = pd.read_csv(path, header=None)
    if df.shape[1] < 2:
        raise ValueError("network_diamond.txt debe tener al menos dos columnas (u,v)")
    G = nx.from_pandas_edgelist(df, source=0, target=1, create_using=nx.Graph())
    return G


def load_guild_network(path: str) -> nx.Graph:
    # Espacios: u w v  (peso intermedio)
    rows = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            u, w, v = parts[0], parts[1], parts[2]
            try:
                w = float(w)
            except Exception:
                w = 1.0
            rows.append((u, v, w))
    df = pd.DataFrame(rows, columns=["u", "v", "weight"])
    G = nx.from_pandas_edgelist(df, source="u", target="v", edge_attr="weight", create_using=nx.Graph())
    return G


def load_network(path: str, fmt: str) -> nx.Graph:
    fmt = fmt.lower()
    if fmt == "string":
        return load_string_network(path)
    elif fmt == "diamond":
        return load_diamond_network(path)
    elif fmt == "guild":
        return load_guild_network(path)
    else:
        raise ValueError("--network-format debe ser uno de: string | diamond | guild")


# ------------------ DIAMOnD ------------------

def diamond_rank(G: nx.Graph, seeds: Set[str], k: int = 200) -> pd.DataFrame:
    """
    Implementación didáctica de DIAMOnD:
      1) Módulo inicial M = semillas
      2) Iterativamente, para cada nodo no en M, calcula p = P(X>=k) con hipergeométrica:
           - N = |V| - 1
           - K = |M|
           - n = grado(v)
           - k = vecinos_en_M(v)
         Selecciona el nodo con p menor y lo añade a M.
      3) Repite hasta añadir k nodos o quedarte sin candidatos.
    Devuelve un DataFrame con el orden de incorporación y p-valor.
    """
    V = set(G.nodes())
    M = [n for n in seeds if n in V]  # conserva orden si quieres
    Mset = set(M)
    if len(Mset) == 0:
        raise ValueError("Ninguna semilla está en la red. ¿Estás usando la red adecuada para HUGO?")

    ranking = []
    step = 0
    while step < k:
        candidates = list(V - Mset)
        if not candidates:
            break

        N = len(V) - 1
        K = len(Mset)
        best_node, best_p, best_k = None, 1.0, 0

        for v in candidates:
            neigh = set(G.neighbors(v))
            d = len(neigh)
            if d == 0:
                continue
            kv = len(neigh & Mset)
            # p-valor hipergeométrico (cola)
            p = hypergeom_sf(kv, N, K, d)
            if p < best_p:
                best_node, best_p, best_k = v, p, kv

        if best_node is None:
            break

        Mset.add(best_node)
        step += 1
        ranking.append({
            "rank": step,
            "node": best_node,
            "k_links_to_module": best_k,
            "p_value": best_p
        })

    df = pd.DataFrame(ranking)
    return df


def read_seeds(path: str) -> List[str]:
    with open(path, "r", encoding="utf-8") as f:
        seeds = [l.strip() for l in f if l.strip()]
    return seeds


def plot_top_bar(df: pd.DataFrame, out_png: str, top: int = 20):
    import matplotlib.pyplot as plt

    d = df.head(top).copy()
    d["score"] = -np.log10(d["p_value"].clip(lower=1e-300))
    d = d.sort_values("score", ascending=True)

    plt.figure(figsize=(10, 8))
    plt.barh(d["node"], d["score"])
    plt.xlabel("-log10(p)")
    plt.ylabel("Nodo")
    plt.title("Top nodos DIAMOnD")
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()


def save_subgraph(G: nx.Graph, seeds: Set[str], top_nodes: Iterable[str], out_edge: str):
    keep = set(seeds) | set(top_nodes)
    S = G.subgraph(keep).copy()
    nx.write_edgelist(S, out_edge, data=["weight"])


def main():
    ap = argparse.ArgumentParser(description="Propagación en redes con DIAMOnD (docente).")
    ap.add_argument("--network", required=True, help="Ruta de la red")
    ap.add_argument("--network-format", required=True, choices=["string", "diamond", "guild"],
                    help="Formato de red: string | diamond | guild")
    ap.add_argument("--seeds", required=True, help="Fichero con semillas (una por línea)")
    ap.add_argument("--k", type=int, default=200, help="Nº de nodos a añadir con DIAMOnD")
    ap.add_argument("--outdir", default="results", help="Directorio de salida")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d-%H%M%S")

    # Cargar red
    G = load_network(args.network, args.network_format)
    # Cargar semillas
    seeds = set(read_seeds(args.seeds))

    # Filtrar semillas existentes en red
    seeds_in = [s for s in seeds if s in G]
    seeds_out = list(set(seeds) - set(seeds_in))

    # Ejecutar DIAMOnD
    df_rank = diamond_rank(G, set(seeds_in), k=args.k)

    # Guardar ranking
    out_tsv = os.path.join(args.outdir, "diamond_ranking.tsv")
    df_rank.to_csv(out_tsv, sep="\t", index=False)

    # Guardar metadatos
    meta = {
        "timestamp": ts,
        "network": args.network,
        "network_format": args.network_format,
        "nodes": int(G.number_of_nodes()),
        "edges": int(G.number_of_edges()),
        "seeds_input": sorted(list(seeds)),
        "seeds_in_network": sorted(list(seeds_in)),
        "seeds_not_found": sorted(list(seeds_out)),
        "k": args.k,
        "scipy_available": _HAVE_SCIPY
    }
    out_json = os.path.join(args.outdir, "diamond_params.json")
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2)

    # Figura top 20
    out_png = os.path.join(args.outdir, "top20_barplot.png")
    if not df_rank.empty:
        plot_top_bar(df_rank, out_png, top=20)

    # Subred top 50
    top50 = df_rank.head(50)["node"].tolist()
    out_edges = os.path.join(args.outdir, "subgraph_top50.edgelist")
    save_subgraph(G, set(seeds_in), top50, out_edges)

    print(f"[OK] DIAMOnD completado.")
    print(f"Semillas en red: {len(seeds_in)} / {len(seeds)}  —  ausentes: {len(seeds_out)}")
    print(f"Ranking:   {out_tsv}")
    print(f"Metadatos: {out_json}")
    if not df_rank.empty:
        print(f"Figura:    {out_png}")
        print(f"Subred:    {out_edges}")
    else:
        print("[Aviso] No se generó ranking (¿semillas no conectan o red vacía?).")


if __name__ == "__main__":
    main()
