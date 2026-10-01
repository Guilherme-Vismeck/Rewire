#!/usr/bin/env python3
"""Exporta as redes em GraphML para o Cytoscape (EXPORT_NETWORK)."""
import argparse
from pathlib import Path

import networkx as nx
import pandas as pd


def group_graph(path):
    df = pd.read_csv(path, sep="\t")
    g = nx.Graph()
    for row in df.itertuples(index=False):
        g.add_edge(
            str(row.gene_a), str(row.gene_b),
            weight=float(row.r), abs_weight=abs(float(row.r)),
            sign=str(row.sign), fdr=float(row.fdr),
        )
    return g


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--network_a", required=True)
    p.add_argument("--network_b", required=True)
    p.add_argument("--edges", required=True, help="differential_edges.tsv")
    p.add_argument("--genes", required=True, help="gene_rewiring.tsv")
    p.add_argument("--label_a", required=True)
    p.add_argument("--label_b", required=True)
    p.add_argument("--outdir", default=".")
    args = p.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    la, lb = args.label_a, args.label_b

    ga, gb = group_graph(args.network_a), group_graph(args.network_b)
    nx.write_graphml(ga, out / f"{la}.graphml")
    nx.write_graphml(gb, out / f"{lb}.graphml")

    edges = pd.read_csv(args.edges, sep="\t")
    genes = pd.read_csv(args.genes, sep="\t", index_col=0)

    g = nx.Graph()
    for row in edges.itertuples(index=False):
        r = row._asdict()
        g.add_edge(
            str(r["gene1"]), str(r["gene2"]),
            status=str(r["status"]),
            **{f"r_{la}": float(r[f"r_{la}"]), f"r_{lb}": float(r[f"r_{lb}"])},
            delta_r=float(r["delta_r"]),
            fdr_diff=float(r["fdr_diff"]),
            dcl=bool(r["dcl"]),
        )

    for node in g.nodes:
        da = ga.degree(node) if node in ga else 0
        db = gb.degree(node) if node in gb else 0
        g.nodes[node]["gene"] = node
        g.nodes[node][f"degree_{la}"] = int(da)
        g.nodes[node][f"degree_{lb}"] = int(db)
        g.nodes[node]["delta_degree"] = int(db - da)
        for col in ["edges_gained", "edges_lost", "sign_flips", "dcl_total"]:
            g.nodes[node][col] = int(genes.loc[node, col]) if node in genes.index else 0

    nx.write_graphml(g, out / "rewiring.graphml")

    print(f"{la}.graphml: {ga.number_of_nodes()} nós, {ga.number_of_edges()} arestas")
    print(f"{lb}.graphml: {gb.number_of_nodes()} nós, {gb.number_of_edges()} arestas")
    print(f"rewiring.graphml: {g.number_of_nodes()} nós, {g.number_of_edges()} arestas")


if __name__ == "__main__":
    main()