#!/usr/bin/env python3
"""Confere módulos e mudança de módulo no dataset sintético."""
import argparse
import sys

import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--modules", required=True)
args = p.parse_args()

df = pd.read_csv(args.modules, sep="\t").set_index("gene")

expected = {
    "GENE1": "changed",
    "GENE2": "lost_module",
    "GENE3": "gained_module",
    "GENE4": "same", "GENE5": "same", "GENE6": "same", "GENE7": "same",
    "GENE8": "no_module",
}

ok = True
for gene, exp in expected.items():
    got = df.loc[gene, "module_status"]
    good = got == exp
    ok &= good
    print(f"{'PASS' if good else 'FAIL'}  {gene}: esperado {exp}, obtido {got}")

bg = df[df.index.str.startswith("BG")]
n_mod = int((bg["module_status"] != "no_module").sum())
good = n_mod == 0
ok &= good
print(f"{'PASS' if good else 'FAIL'}  genes de fundo em algum módulo: {n_mod}")

sys.exit(0 if ok else 1)
