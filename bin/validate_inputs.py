#!/usr/bin/env python3
"""Valida expression.tsv e metadata.tsv antes da análise (VALIDATE_INPUTS)."""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--expression", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--group_col", default="condition")
    p.add_argument("--group_a", required=True)
    p.add_argument("--group_b", required=True)
    p.add_argument("--gene_column", default=None, help="padrão: primeira coluna")
    p.add_argument("--sample_column", default="sample")
    p.add_argument("--report", default="qc_report.tsv")
    args = p.parse_args()

    errors, warnings, info = [], [], []

    # 1. Arquivos existem
    for f in (args.expression, args.metadata):
        if not Path(f).is_file():
            sys.exit(f"ERROR: arquivo não encontrado: {f}")

    # 2. Leitura
    expr = pd.read_csv(args.expression, sep="\t", dtype=str)
    meta = pd.read_csv(args.metadata, sep="\t", dtype=str)

    gene_col = args.gene_column or expr.columns[0]
    if gene_col not in expr.columns:
        sys.exit(f"ERROR: coluna de genes '{gene_col}' não existe em {args.expression}")

    # Amostras duplicadas no cabeçalho (o pandas renomearia silenciosamente)
    with open(args.expression) as fh:
        header = fh.readline().rstrip("\n").split("\t")
    dup_cols = sorted({c for c in header if header.count(c) > 1})
    if dup_cols:
        errors.append(f"nomes de amostra duplicados na matriz: {dup_cols}")

    # 3. Genes: existência e duplicatas
    if expr.shape[0] == 0:
        errors.append("nenhum gene encontrado na matriz")
    dup_genes = expr.loc[expr[gene_col].duplicated(keep=False), gene_col].unique()
    if len(dup_genes) > 0:
        errors.append(f"{len(dup_genes)} genes duplicados, ex.: {list(dup_genes[:5])}")

    matrix_samples = [c for c in expr.columns if c != gene_col]
    if len(matrix_samples) == 0:
        errors.append("nenhuma amostra encontrada na matriz")

    # 4. Valores numéricos, ausentes e infinitos
    values = expr[matrix_samples].apply(pd.to_numeric, errors="coerce")
    non_numeric = int((values.isna() & expr[matrix_samples].notna()).sum().sum())
    if non_numeric > 0:
        errors.append(f"{non_numeric} valores não numéricos na matriz")
    n_missing = int(values.isna().sum().sum())
    if n_missing > 0:
        warnings.append(f"{n_missing} valores ausentes (NA) na matriz")
    n_inf = int(np.isinf(values.to_numpy(dtype=float)).sum())
    if n_inf > 0:
        errors.append(f"{n_inf} valores infinitos na matriz")

    # 5. Variância zero
    zero_var = int((values.var(axis=1, skipna=True).fillna(0) == 0).sum())
    if zero_var > 0:
        warnings.append(f"{zero_var} genes com variância zero (serão removidos no PREPROCESS)")

    # 6. Metadata
    for col in (args.sample_column, args.group_col):
        if col not in meta.columns:
            errors.append(f"coluna '{col}' não existe em {args.metadata}")
    if errors and (args.sample_column not in meta.columns or args.group_col not in meta.columns):
        write_report(args.report, errors, warnings, info)
        sys.exit("ERROR:\n  - " + "\n  - ".join(errors))

    if meta[args.sample_column].duplicated().any():
        errors.append("amostras duplicadas no metadata")

    # 7. Amostras da matriz x metadata
    meta_samples = set(meta[args.sample_column])
    missing_in_meta = [s for s in matrix_samples if s not in meta_samples]
    if missing_in_meta:
        errors.append(
            f"{len(missing_in_meta)} amostras da matriz não estão no metadata: {missing_in_meta[:10]}"
        )
    extra_in_meta = [s for s in meta_samples if s not in set(matrix_samples)]
    if extra_in_meta:
        warnings.append(
            f"{len(extra_in_meta)} amostras do metadata não estão na matriz (serão ignoradas): {extra_in_meta[:10]}"
        )

    # 8. Grupos solicitados e tamanho amostral
    usable = meta[meta[args.sample_column].isin(matrix_samples)]
    counts = usable[args.group_col].value_counts()
    for g in (args.group_a, args.group_b):
        if g not in counts.index:
            errors.append(f"grupo '{g}' não encontrado na coluna '{args.group_col}'")
            continue
        n = int(counts[g])
        info.append(f"grupo '{g}': {n} amostras")
        if n < 15:
            warnings.append(f"grupo '{g}' com n={n} (<15): estimativas de rede altamente instáveis")
        elif n < 30:
            warnings.append(f"grupo '{g}' com n={n} (15-29): análise exploratória")

    info.append(f"{expr.shape[0]} genes x {len(matrix_samples)} amostras")

    write_report(args.report, errors, warnings, info)

    for m in info:
        print(f"INFO: {m}")
    for m in warnings:
        print(f"WARNING: {m}")
    if errors:
        print("\nERROR:")
        for m in errors:
            print(f"  - {m}")
        sys.exit(1)
    print("Validação concluída sem erros.")


def write_report(path, errors, warnings, info):
    rows = (
        [("ERROR", m) for m in errors]
        + [("WARNING", m) for m in warnings]
        + [("INFO", m) for m in info]
    )
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=["level", "message"]).to_csv(path, sep="\t", index=False)


if __name__ == "__main__":
    main()