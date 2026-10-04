#!/usr/bin/env python3
"""Bootstrap da estabilidade das arestas e rewiring estável (BOOTSTRAP)."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

CHANGED = ["GAINED", "LOST", "SIGN_FLIPPED", "STRENGTHENED", "WEAKENED"]
GENE_COLS = ["raw_rewired_edges", "stable_rewired_edges", "mean_support", "stable_fraction"]


def pair_corr(X, gi, gj, idx, spearman):
    """Correlação dos pares (gi, gj) numa reamostragem das amostras (colunas idx)."""
    Xb = X[:, idx]
    if spearman:
        Xb = stats.rankdata(Xb, axis=1)
    Xb = Xb - Xb.mean(axis=1, keepdims=True)
    norm = np.sqrt((Xb ** 2).sum(axis=1))
    with np.errstate(divide="ignore", invalid="ignore"):
        Z = Xb / norm[:, None]
        return (Z[gi] * Z[gj]).sum(axis=1)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--edges", required=True, help="differential_edges.tsv")
    p.add_argument("--expr_a", required=True)
    p.add_argument("--expr_b", required=True)
    p.add_argument("--label_a", required=True)
    p.add_argument("--label_b", required=True)
    p.add_argument("--method", choices=["pearson", "spearman"], default="pearson")
    p.add_argument("--min_cor", type=float, default=0.6)
    p.add_argument("--bootstrap", type=int, default=100)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--stability_threshold", type=float, default=0.8)
    p.add_argument("--outdir", default=".")
    args = p.parse_args()

    la, lb = args.label_a, args.label_b
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    edges = pd.read_csv(args.edges, sep="\t")
    stab_cols = [f"stability_{la}", f"stability_{lb}"]
    edge_cols = ["gene1", "gene2", "status"] + stab_cols + [
        "bootstrap_support", "raw_rewiring", "stable_rewiring"]
    stable_cols = ["gene1", "gene2", "status", "delta_r", "fdr_diff"] + stab_cols + [
        "bootstrap_support"]

    if edges.empty:
        pd.DataFrame(columns=edge_cols).to_csv(out / "edge_stability.tsv", sep="\t", index=False)
        pd.DataFrame(columns=stable_cols).to_csv(out / "stable_rewiring.tsv", sep="\t", index=False)
        pd.DataFrame(columns=["gene"] + GENE_COLS).to_csv(
            out / "gene_stability.tsv", sep="\t", index=False)
        print("Nenhuma aresta candidata; arquivos vazios gerados.")
        return

    genes = sorted(set(edges["gene1"].astype(str)) | set(edges["gene2"].astype(str)))
    pos = {g: i for i, g in enumerate(genes)}
    gi = edges["gene1"].astype(str).map(pos).to_numpy()
    gj = edges["gene2"].astype(str).map(pos).to_numpy()

    rng = np.random.default_rng(args.seed)
    spearman = args.method == "spearman"

    for label, path in [(la, args.expr_a), (lb, args.expr_b)]:
        expr = pd.read_csv(path, sep="\t", index_col=0)
        expr.index = expr.index.astype(str)
        X = expr.loc[genes].to_numpy(dtype=float)
        n = X.shape[1]
        count = np.zeros(len(edges), dtype=int)
        for _ in range(args.bootstrap):
            idx = rng.integers(0, n, size=n)
            r = pair_corr(X, gi, gj, idx, spearman)
            count += (np.abs(r) >= args.min_cor)
        edges[f"stability_{label}"] = count / args.bootstrap

    sa, sb = edges[f"stability_{la}"], edges[f"stability_{lb}"]
    status = edges["status"]
    edges["bootstrap_support"] = np.where(
        status == "LOST", np.minimum(sa, 1 - sb),
        np.where(status == "GAINED", np.minimum(sb, 1 - sa), np.minimum(sa, sb)))

    edges["raw_rewiring"] = edges["dcl"].astype(bool) & status.isin(CHANGED)
    edges["stable_rewiring"] = edges["raw_rewiring"] & (
        edges["bootstrap_support"] >= args.stability_threshold)

    edges.sort_values("bootstrap_support", ascending=False)[edge_cols].to_csv(
        out / "edge_stability.tsv", sep="\t", index=False, float_format="%.4g")
    (edges[edges["stable_rewiring"]]
        .sort_values(["bootstrap_support", "fdr_diff"], ascending=[False, True])[stable_cols]
        .to_csv(out / "stable_rewiring.tsv", sep="\t", index=False, float_format="%.4g"))

    raw = edges[edges["raw_rewiring"]]
    long = pd.concat([
        raw[["gene1", "stable_rewiring", "bootstrap_support"]].rename(columns={"gene1": "gene"}),
        raw[["gene2", "stable_rewiring", "bootstrap_support"]].rename(columns={"gene2": "gene"}),
    ])
    if long.empty:
        gene = pd.DataFrame(columns=GENE_COLS)
        gene.index.name = "gene"
    else:
        gene = long.groupby("gene").agg(
            raw_rewired_edges=("stable_rewiring", "size"),
            stable_rewired_edges=("stable_rewiring", "sum"),
            mean_support=("bootstrap_support", "mean"))
        gene["stable_fraction"] = gene["stable_rewired_edges"] / gene["raw_rewired_edges"]
        gene = gene.sort_values(["stable_rewired_edges", "mean_support"], ascending=False)
    gene.to_csv(out / "gene_stability.tsv", sep="\t", float_format="%.4g")

    print(f"Repetições de bootstrap: {args.bootstrap} (seed {args.seed})")
    print(f"Arestas candidatas avaliadas: {len(edges)}")
    print(f"Arestas com rewiring bruto (DCL): {int(edges['raw_rewiring'].sum())}")
    print(f"Arestas com rewiring estável (suporte >= {args.stability_threshold}): "
          f"{int(edges['stable_rewiring'].sum())}")


if __name__ == "__main__":
    main()
