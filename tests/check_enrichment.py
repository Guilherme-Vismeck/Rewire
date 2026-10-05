#!/usr/bin/env python3
"""Confere o enriquecimento com conjuntos próprios no dataset sintético."""
import argparse
import sys

import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--enrichment", required=True)
args = p.parse_args()

df = pd.read_csv(args.enrichment, sep="\t")
ok = True

for gene_list in ["rewired_stable", "RW_only", "gained_connections",
                  "lost_connections", "module_changed"]:
    row = df[(df["gene_list"] == gene_list) & (df["term"] == "REWIRED_PROGRAM")]
    good = len(row) == 1 and row.iloc[0]["fdr"] < 0.05
    ok &= good
    detail = f"sobreposição {int(row.iloc[0]['overlap'])}, FDR {row.iloc[0]['fdr']:.2g}" if len(row) else "ausente"
    print(f"{'PASS' if good else 'FAIL'}  {gene_list}: REWIRED_PROGRAM enriquecido ({detail})")

row = df[(df["gene_list"] == "RW_only") & (df["term"] == "REWIRED_PROGRAM")]
good = len(row) == 1 and int(row.iloc[0]["overlap"]) == 5 and abs(row.iloc[0]["fold_enrichment"] - 108 / 5) < 0.1
ok &= good
print(f"{'PASS' if good else 'FAIL'}  RW_only: 5 de 5 genes e fold enrichment ~ 21,6")

wrong = df[df["term"].isin(["DECOY_BACKGROUND", "STABLE_PAIR"])]
good = len(wrong) == 0
ok &= good
print(f"{'PASS' if good else 'FAIL'}  controles negativos enriquecidos indevidamente: {len(wrong)}")

sys.exit(0 if ok else 1)
