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

bg = df[df.index.str.startswith("BG")]
bad = bg[bg["category"] != "unchanged"]
good = len(bad) == 0
ok &= good
print(f"{'PASS' if good else 'FAIL'}  genes de fundo fora de 'unchanged': {len(bad)}")

sys.exit(0 if ok else 1)
