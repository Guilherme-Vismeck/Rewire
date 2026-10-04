#!/usr/bin/env python3
"""Cruza expressão diferencial com rewiring estável (INTEGRATE_DE_RW)."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--de", required=True, help="DE_results.tsv")
    p.add_argument("--gene_stability", required=True, help="gene_stability.tsv")
    p.add_argument("--min_rewired_edges", type=int, default=1,
                   help="arestas de rewiring estável para considerar o gene rewired")
    p.add_argument("--outdir", default=".")
    args = p.parse_args()

    de = pd.read_csv(args.de, sep="\t")
    gs = pd.read_csv(args.gene_stability, sep="\t")
    df = de.merge(gs, on="gene", how="left")

    for col in ["raw_rewired_edges", "stable_rewired_edges"]:
        df[col] = df[col].fillna(0).astype(int)
    df["DEG"] = df["DEG"].fillna(False).astype(bool)
    df["rewired"] = df["stable_rewired_edges"] >= args.min_rewired_edges

    df["category"] = np.select(
        [df["DEG"] & df["rewired"], ~df["DEG"] & df["rewired"], df["DEG"] & ~df["rewired"]],
        ["DEG+RW", "RW_only", "DEG_only"],
        default="unchanged",
    )

    order = {"RW_only": 0, "DEG+RW": 1, "DEG_only": 2, "unchanged": 3}
    df = (df.assign(_o=df["category"].map(order))
            .sort_values(["_o", "stable_rewired_edges", "pvalue"], ascending=[True, False, True])
            .drop(columns="_o"))

    summary = (df["category"].value_counts()
               .reindex(["DEG+RW", "RW_only", "DEG_only", "unchanged"], fill_value=0)
               .rename_axis("category").reset_index(name="n_genes"))

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "de_vs_rewiring.tsv", sep="\t", index=False, float_format="%.4g")
    summary.to_csv(out / "category_summary.tsv", sep="\t", index=False)

    print("Genes por categoria:")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
