#!/usr/bin/env python3
"""Confere o rewiring estável no dataset sintético."""
import argparse
import sys

import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--stable", required=True)
p.add_argument("--edge_stability", required=True)
args = p.parse_args()

stable = pd.read_csv(args.stable, sep="\t")
found = {frozenset((a, b)): s for a, b, s in zip(stable["gene1"], stable["gene2"], stable["status"])}
expected = {
    frozenset(("GENE1", "GENE2")): "LOST",
    frozenset(("GENE1", "GENE3")): "GAINED",
    frozenset(("GENE4", "GENE5")): "SIGN_FLIPPED",
}

ok = True
for pair, status in expected.items():
    got = found.get(pair, "AUSENTE")
    good = got == status
    ok &= good
    print(f"{'PASS' if good else 'FAIL'}  {'-'.join(sorted(pair))}: esperado {status} estável, obtido {got}")

extras = set(found) - set(expected)
good = len(extras) == 0
ok &= good
print(f"{'PASS' if good else 'FAIL'}  arestas estáveis extras: {len(extras)}")

es = pd.read_csv(args.edge_stability, sep="\t")
mask = es.apply(lambda r: frozenset((r["gene1"], r["gene2"])) == frozenset(("GENE6", "GENE7")), axis=1)
row = es[mask].iloc[0]
good = row["bootstrap_support"] >= 0.9 and not bool(row["stable_rewiring"])
ok &= good
print(f"{'PASS' if good else 'FAIL'}  GENE6-GENE7: suporte {row['bootstrap_support']} e não marcada como rewiring")

sys.exit(0 if ok else 1)
