#!/usr/bin/env python3
"""Métricas topológicas e neighbor turnover entre dois grupos (TOPOLOGY)."""
import argparse
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd

METRICS = ["degree", "strength", "betweenness", "closeness", "eigenvector", "clustering"]


def build_graph(path, genes):
    g = nx.Graph()
    g.add_nodes_from(genes)  # genes sem arestas entram como nós isolados
    df = pd.read_csv(path, sep="\t")
    for row in df.itertuples(index=False):
        w = abs(float(row.r))
        g.add_edge(str(row.gene_a), str(row.gene_b), abs_weight=w, distance=1.0 / w)
    return g


def compute_metrics(g):
    try:
        eig = {n: abs(v) for n, v in
               nx.eigenvector_centrality_numpy(g, weight="abs_weight").items()}
    except Exception:
        eig = {n: 0.0 for n in g}
    return pd.DataFrame({
        "degree": dict(g.degree()),
        "strength": dict(g.degree(weight="abs_weight")),
        "betweenness": nx.betweenness_centrality(g, weight="distance"),
        "closeness": nx.closeness_centrality(g, distance="distance"),
        "eigenvector": eig,
        "clustering": nx.clustering(g, weight="abs_weight"),
    })


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--network_a", required=True)
    p.add_argument("--network_b", required=True)
    p.add_argument("--expression", required=True, help="filtered_expression.tsv (lista de genes)")
    p.add_argument("--label_a", required=True)
    p.add_argument("--label_b", required=True)
    p.add_argument("--outdir", default=".")
    args = p.parse_args()

    la, lb = args.label_a, args.label_b
    genes = pd.read_csv(args.expression, sep="\t", usecols=[0]).iloc[:, 0].astype(str).tolist()

    ga = build_graph(args.network_a, genes)
    gb = build_graph(args.network_b, genes)

    ma, mb = compute_metrics(ga), compute_metrics(gb)
    table = ma.add_suffix(f"_{la}").join(mb.add_suffix(f"_{lb}"))
    for m in METRICS:
        table[f"delta_{m}"] = table[f"{m}_{lb}"] - table[f"{m}_{la}"]
    table.index.name = "gene"

    rows = []
    for gene in genes:
        na, nb = set(ga[gene]), set(gb[gene])
        union, inter = na | nb, na & nb
        jac = len(inter) / len(union) if union else np.nan
        rows.append({
            "gene": gene,
            f"neighbors_{la}": len(na),
            f"neighbors_{lb}": len(nb),
            "neighbors_shared": len(inter),
            "neighbors_gained": len(nb - na),
            "neighbors_lost": len(na - nb),
            "neighbor_jaccard": jac,
            "neighbor_turnover": 1 - jac if union else np.nan,
        })
    turnover = (pd.DataFrame(rows).set_index("gene")
                .sort_values("neighbor_turnover", ascending=False))

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    table.sort_values(f"strength_{lb}", ascending=False).to_csv(
        out / "network_metrics.tsv", sep="\t", float_format="%.5g")
    turnover.to_csv(out / "neighbor_turnover.tsv", sep="\t", float_format="%.5g")

    n_conn = int(turnover["neighbor_turnover"].notna().sum())
    print(f"Genes analisados: {len(genes)}")
    print(f"Genes com vizinhos em pelo menos um grupo: {n_conn}")
    print(f"Genes com turnover = 1 (vizinhos totalmente diferentes): "
          f"{int((turnover['neighbor_turnover'] == 1).sum())}")


if __name__ == "__main__":
    main()
