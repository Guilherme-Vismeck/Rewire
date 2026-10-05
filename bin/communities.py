#!/usr/bin/env python3
"""Módulos (Louvain) em cada rede e mudança de módulo por gene (COMMUNITIES)."""
import argparse
from pathlib import Path

import networkx as nx
import pandas as pd
from networkx.algorithms.community import louvain_communities

STATUS_ORDER = ["changed", "lost_module", "gained_module", "same", "no_module"]


def build_graph(path, genes):
    g = nx.Graph()
    g.add_nodes_from(genes)
    df = pd.read_csv(path, sep="\t")
    for row in df.itertuples(index=False):
        g.add_edge(str(row.gene_a), str(row.gene_b), abs_weight=abs(float(row.r)))
    return g


def detect_modules(g, resolution, seed, min_size):
    comms = louvain_communities(g, weight="abs_weight", resolution=resolution, seed=seed)
    comms = [set(c) for c in comms if len(c) >= min_size]
    comms.sort(key=lambda c: (-len(c), min(c)))  # ordem determinística
    return {f"M{i + 1}": c for i, c in enumerate(comms)}


def membership(modules):
    return {gene: (name, members) for name, members in modules.items() for gene in members}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--network_a", required=True)
    p.add_argument("--network_b", required=True)
    p.add_argument("--expression", required=True, help="filtered_expression.tsv (lista de genes)")
    p.add_argument("--label_a", required=True)
    p.add_argument("--label_b", required=True)
    p.add_argument("--resolution", type=float, default=1.0)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--min_module_size", type=int, default=2)
    p.add_argument("--module_overlap_threshold", type=float, default=0.5)
    p.add_argument("--outdir", default=".")
    args = p.parse_args()

    la, lb = args.label_a, args.label_b
    genes = pd.read_csv(args.expression, sep="\t", usecols=[0]).iloc[:, 0].astype(str).tolist()

    ga, gb = build_graph(args.network_a, genes), build_graph(args.network_b, genes)
    mods_a = detect_modules(ga, args.resolution, args.seed, args.min_module_size)
    mods_b = detect_modules(gb, args.resolution, args.seed, args.min_module_size)
    mem_a, mem_b = membership(mods_a), membership(mods_b)

    rows = []
    for gene in genes:
        ia, ib = mem_a.get(gene), mem_b.get(gene)
        jac = None
        if ia is None and ib is None:
            status = "no_module"
        elif ia is None:
            status = "gained_module"
        elif ib is None:
            status = "lost_module"
        else:
            jac = len(ia[1] & ib[1]) / len(ia[1] | ib[1])
            status = "same" if jac >= args.module_overlap_threshold else "changed"
        rows.append({
            "gene": gene,
            f"module_{la}": ia[0] if ia else None,
            f"module_{lb}": ib[0] if ib else None,
            f"module_size_{la}": len(ia[1]) if ia else 0,
            f"module_size_{lb}": len(ib[1]) if ib else 0,
            "module_jaccard": jac,
            "module_status": status,
        })

    table = pd.DataFrame(rows)
    table["_o"] = table["module_status"].map({s: i for i, s in enumerate(STATUS_ORDER)})
    table = table.sort_values(["_o", "gene"]).drop(columns="_o")

    summary = pd.DataFrame(
        [(label, name, len(members), ",".join(sorted(members)))
         for label, mods in [(la, mods_a), (lb, mods_b)] for name, members in mods.items()],
        columns=["group", "module", "n_genes", "genes"],
    )

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    table.to_csv(out / "gene_modules.tsv", sep="\t", index=False, float_format="%.4g")
    summary.to_csv(out / "module_summary.tsv", sep="\t", index=False)

    print(f"[{la}] {len(mods_a)} módulos (tamanho >= {args.min_module_size})")
    print(f"[{lb}] {len(mods_b)} módulos (tamanho >= {args.min_module_size})")
    print("Genes por status:")
    print(table["module_status"].value_counts().reindex(STATUS_ORDER, fill_value=0).to_string())


if __name__ == "__main__":
    main()
