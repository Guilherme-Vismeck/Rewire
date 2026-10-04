#!/usr/bin/env python3
"""Confere as categorias DE x rewiring no dataset sintético."""
import argparse
import sys

import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--table", required=True)
args = p.parse_args()

df = pd.read_csv(args.table, sep="\t").set_index("gene")

expected = {
    "GENE1": "RW_only", "GENE2": "RW_only", "GENE3": "RW_only",
    "GENE4": "RW_only", "GENE5": "RW_only",
    "GENE6": "unchanged", "GENE7": "unchanged",
    "GENE8": "DEG_only",
}

ok = True
for gene, exp in expected.items():
    got = df.loc[gene, "category"]
    good = got == exp
    ok &= good
    print(f"{'PASS' if good else 'FAIL'}  {gene}: esperado {exp}, obtido {got}")

# Genes de fundo são ruído: nenhum pode ter rewiring estável.
bg = df[df.index.str.startswith("BG")]
n_rw = int(bg["rewired"].sum())
good = n_rw == 0
ok &= good
print(f"{'PASS' if good else 'FAIL'}  genes de fundo com rewiring estável: {n_rw}")

# Falsos positivos de DE são esperados em pequena fração (limma com FDR de 5%).
max_false_deg = max(1, int(0.05 * len(bg)))
n_deg = int(bg["DEG"].sum())
good = n_deg <= max_false_deg
ok &= good
print(f"{'PASS' if good else 'FAIL'}  genes de fundo marcados como DEG: {n_deg} "
      f"(tolerância de ruído: até {max_false_deg} de {len(bg)})")
if n_deg:
    print("      ruído:", ", ".join(bg.index[bg["DEG"]]))

sys.exit(0 if ok else 1)
