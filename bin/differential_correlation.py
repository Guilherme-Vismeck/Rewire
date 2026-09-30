#!/usr/bin/env python3
"""Correlação diferencial entre dois grupos e classificação das arestas (DIFFERENTIAL_NETWORK)."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

CHANGE_COLS = ["GAINED", "LOST", "SIGN_FLIPPED", "STRENGTHENED", "WEAKENED"]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--corr_a", required=True, help="*_correlations.tsv.gz do grupo A")
    p.add_argument("--corr_b", required=True, help="*_correlations.tsv.gz do grupo B")
    p.add_argument("--label_a", required=True)
    p.add_argument("--label_b", required=True)
    p.add_argument("--n_a", type=int, required=True, help="nº de amostras do grupo A")
    p.add_argument("--n_b", type=int, required=True, help="nº de amostras do grupo B")
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--outdir", default=".")
    args = p.parse_args()

    cols = ["gene_a", "gene_b", "r", "present"]
    a = pd.read_csv(args.corr_a, sep="\t", usecols=cols)
    b = pd.read_csv(args.corr_b, sep="\t", usecols=cols)
    df = a.merge(b, on=["gene_a", "gene_b"], suffixes=("_a", "_b"))
    if not (len(df) == len(a) == len(b)):
        raise SystemExit("ERROR: as duas redes não cobrem o mesmo conjunto de pares de genes")
    # gene_a/gene_b aqui são os dois genes do par, não os grupos
    df = df.rename(columns={"gene_a": "gene1", "gene_b": "gene2"})

    # Teste de Fisher: diferença entre duas correlações independentes
    clip = 0.999999
    za = np.arctanh(np.clip(df["r_a"], -clip, clip))
    zb = np.arctanh(np.clip(df["r_b"], -clip, clip))
    se = np.sqrt(1 / (args.n_a - 3) + 1 / (args.n_b - 3))
    df["delta_r"] = df["r_b"] - df["r_a"]
    df["z_diff"] = (zb - za) / se
    df["p_diff"] = 2 * stats.norm.sf(np.abs(df["z_diff"]))
    df["fdr_diff"] = multipletests(df["p_diff"], method="fdr_bh")[1]
    df["dcl"] = df["fdr_diff"] < args.fdr

    # Classificação
    pa, pb = df["present_a"], df["present_b"]
    flipped = np.sign(df["r_a"]) != np.sign(df["r_b"])
    stronger = df["r_b"].abs() > df["r_a"].abs()
    conditions = [
        pa & pb & flipped,
        pa & ~pb,
        ~pa & pb,
        pa & pb & df["dcl"] & stronger,
        pa & pb & df["dcl"] & ~stronger,
        pa & pb,
    ]
    choices = ["SIGN_FLIPPED", "LOST", "GAINED", "STRENGTHENED", "WEAKENED", "PRESERVED"]
    df["status"] = np.select(conditions, choices, default="NOT_PRESENT")

    la, lb = args.label_a, args.label_b
    df = df.rename(columns={
        "r_a": f"r_{la}", "r_b": f"r_{lb}",
        "present_a": f"present_{la}", "present_b": f"present_{lb}",
    })

    edges = df[df["status"] != "NOT_PRESENT"].sort_values("fdr_diff")
    sig = edges[edges["dcl"]]

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    edges.to_csv(out / "differential_edges.tsv", sep="\t", index=False)
    sig[sig["status"] == "GAINED"].to_csv(out / "gained_edges.tsv", sep="\t", index=False)
    sig[sig["status"] == "LOST"].to_csv(out / "lost_edges.tsv", sep="\t", index=False)
    sig[sig["status"] == "SIGN_FLIPPED"].to_csv(out / "sign_flips.tsv", sep="\t", index=False)

    # Resumo por gene (candidatos a DCGs): só arestas com diferença significativa
    if len(sig) > 0:
        long = pd.concat([
            sig[["gene1", "status"]].rename(columns={"gene1": "gene"}),
            sig[["gene2", "status"]].rename(columns={"gene2": "gene"}),
        ])
        genes = (long.groupby(["gene", "status"]).size().unstack(fill_value=0)
                 .reindex(columns=CHANGE_COLS, fill_value=0))
    else:
        genes = pd.DataFrame(columns=CHANGE_COLS)
        genes.index.name = "gene"
    genes.columns = ["edges_gained", "edges_lost", "sign_flips",
                     "edges_strengthened", "edges_weakened"]
    genes["dcl_total"] = genes.sum(axis=1)
    genes = genes.sort_values("dcl_total", ascending=False)
    genes.to_csv(out / "gene_rewiring.tsv", sep="\t")

    print(f"Pares comparados: {len(df)}")
    print("Arestas por status:")
    print(edges["status"].value_counts().to_string())
    print(f"Arestas com diferença significativa (FDR < {args.fdr}): {len(sig)}")
    print(f"Genes com pelo menos uma aresta diferencial: {len(genes)}")


if __name__ == "__main__":
    main()