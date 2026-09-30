#!/usr/bin/env python3
"""Filtra genes, trata NAs, seleciona genes variáveis e separa os grupos (PREPROCESS)."""
import argparse
from pathlib import Path

import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--expression", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--group_col", default="condition")
    p.add_argument("--group_a", required=True)
    p.add_argument("--group_b", required=True)
    p.add_argument("--sample_column", default="sample")
    p.add_argument("--top_variable_genes", type=int, default=3000)
    p.add_argument("--max_missing", type=float, default=0.2,
                   help="fração máxima de NAs por gene")
    p.add_argument("--outdir", default=".")
    args = p.parse_args()

    expr = pd.read_csv(args.expression, sep="\t", index_col=0).apply(pd.to_numeric)
    meta = pd.read_csv(args.metadata, sep="\t", dtype=str)

    groups = {}
    for g in (args.group_a, args.group_b):
        samples = meta.loc[meta[args.group_col] == g, args.sample_column]
        groups[g] = [s for s in samples if s in expr.columns]
    expr = expr[groups[args.group_a] + groups[args.group_b]]
    n_start = expr.shape[0]

    # 1. genes com muitos valores ausentes
    expr = expr[expr.isna().mean(axis=1) <= args.max_missing].copy()
    n_after_na = expr.shape[0]

    # 2. imputação pela mediana do gene dentro de cada grupo
    for cols in groups.values():
        expr[cols] = expr[cols].apply(lambda row: row.fillna(row.median()), axis=1)
    expr = expr.dropna()

    # 3. remove genes com variância zero em qualquer grupo
    var = pd.DataFrame({g: expr[cols].var(axis=1) for g, cols in groups.items()})
    keep = (var > 0).all(axis=1)
    expr, var = expr[keep], var[keep]
    n_after_var = expr.shape[0]

    # 4. seleciona os genes mais variáveis (média das variâncias intra-grupo)
    if expr.shape[0] > args.top_variable_genes:
        top = set(var.mean(axis=1).nlargest(args.top_variable_genes).index)
        expr = expr[expr.index.isin(top)]

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    expr.to_csv(out / "filtered_expression.tsv", sep="\t")
    for g, cols in groups.items():
        expr[cols].to_csv(out / f"{g}_expression.tsv", sep="\t")

    print(f"Genes iniciais: {n_start}")
    print(f"Após filtro de NAs: {n_after_na}")
    print(f"Após filtro de variância zero: {n_after_var}")
    print(f"Genes finais: {expr.shape[0]}")
    for g, cols in groups.items():
        print(f"Grupo '{g}': {len(cols)} amostras")


if __name__ == "__main__":
    main()