#!/usr/bin/env python3
"""Confere topologia e neighbor turnover no dataset sintético."""
import argparse
import sys

import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--turnover", required=True)
p.add_argument("--metrics", required=True)
args = p.parse_args()

turn = pd.read_csv(args.turnover, sep="\t", index_col=0)
met = pd.read_csv(args.metrics, sep="\t", index_col=0)

expected_turnover = {
    "GENE1": 1.0, "GENE2": 1.0, "GENE3": 1.0,
    "GENE4": 0.0, "GENE5": 0.0, "GENE6": 0.0, "GENE7": 0.0,
}
ok = True
for gene, exp in expected_turnover.items():
    got = turn.loc[gene, "neighbor_turnover"]
    good = abs(got - exp) < 1e-9
    ok &= good
    print(f"{'PASS' if good else 'FAIL'}  {gene}: turnover esperado {exp}, obtido {got}")

for gene, col, exp in [("GENE1", "degree_control", 1), ("GENE1", "degree_case", 1),
                       ("GENE2", "degree_case", 0), ("GENE3", "degree_control", 0)]:
    got = met.loc[gene, col]
    good = got == exp
    ok &= good
    print(f"{'PASS' if good else 'FAIL'}  {gene}: {col} esperado {exp}, obtido {got}")

sys.exit(0 if ok else 1)
