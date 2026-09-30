#!/usr/bin/env python3
"""Compara differential_edges.tsv com a verdade conhecida do dataset sintético."""
import argparse
import sys

import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--edges", required=True)
p.add_argument("--truth", required=True)
args = p.parse_args()

edges = pd.read_csv(args.edges, sep="\t")
truth = pd.read_csv(args.truth, sep="\t")

ok = True
expected_pairs = set()
for _, row in truth.iterrows():
    a, b = row["gene_a"], row["gene_b"]
    expected_pairs.add(frozenset((a, b)))
    m = edges[((edges["gene1"] == a) & (edges["gene2"] == b)) |
              ((edges["gene1"] == b) & (edges["gene2"] == a))]
    got = m["status"].iloc[0] if len(m) else "MISSING"
    flag = "PASS" if got == row["expected"] else "FAIL"
    ok &= got == row["expected"]
    print(f"{flag}  {a}-{b}: esperado {row['expected']}, obtido {got}")

extras = edges[~edges.apply(lambda r: frozenset((r["gene1"], r["gene2"])) in expected_pairs, axis=1)]
extras = extras[extras["dcl"]]
print(f"Arestas diferenciais extras (fora da verdade conhecida): {len(extras)}")

sys.exit(0 if ok else 1)