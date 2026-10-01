#!/usr/bin/env python3
"""Constrói a rede de coexpressão de um grupo (BUILD_NETWORK)."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--expression", required=True, help="matriz genes x amostras de um grupo")
    p.add_argument("--label", required=True, help="nome do grupo, usado nos arquivos de saída")
    p.add_argument("--method", choices=["pearson", "spearman"], default="pearson")
    p.add_argument("--min_cor", type=float, default=0.6)
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--outdir", default=".")
    args = p.parse_args()

    expr = pd.read_csv(args.expression, sep="\t", index_col=0)
    genes = expr.index.to_numpy()
    X = expr.to_numpy(dtype=float)
    n = X.shape[1]

    if args.method == "spearman":
        X = np.apply_along_axis(stats.rankdata, 1, X)

    R = np.corrcoef(X)
    iu = np.triu_indices(len(genes), k=1)
    r = R[iu]

    # p-valor bicaudal pela distribuição t com n-2 graus de liberdade
    r_safe = np.clip(r, -0.999999999, 0.999999999)
    t = r_safe * np.sqrt((n - 2) / (1 - r_safe**2))
    pval = 2 * stats.t.sf(np.abs(t), df=n - 2)
    fdr = multipletests(pval, method="fdr_bh")[1]

    df = pd.DataFrame({
        "gene_a": genes[iu[0]],
        "gene_b": genes[iu[1]],
        "r": r,
        "pvalue": pval,
        "fdr": fdr,
    })
    df["present"] = (df["r"].abs() >= args.min_cor) & (df["fdr"] < args.fdr)
    df["sign"] = np.where(df["r"] >= 0, "positive", "negative")

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / f"{args.label}_correlations.tsv.gz", sep="\t", index=False)
    df[df["present"]].to_csv(out / f"{args.label}_network.tsv", sep="\t", index=False)

    print(f"[{args.label}] {len(genes)} genes, {n} amostras, método {args.method}")
    print(f"[{args.label}] {len(df)} pares testados, {int(df['present'].sum())} arestas "
          f"(|r| >= {args.min_cor}, FDR < {args.fdr})")


if __name__ == "__main__":
    main()