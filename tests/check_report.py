#!/usr/bin/env python3
"""Confere o relatório HTML gerado para o dataset sintético."""
import argparse
import sys
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument("--report", required=True)
args = p.parse_args()

doc = Path(args.report).read_text(encoding="utf-8")
n_figs = doc.count("data:image/png;base64")

checks = [
    ("título presente", "RewireNF report" in doc),
    ("tamanho > 30 kB", len(doc) > 30_000),
    (f"pelo menos 7 figuras embutidas (encontradas: {n_figs})", n_figs >= 7),
    ("genes RW_only listados", "RW_only" in doc and "GENE1" in doc),
    ("enriquecimento presente", "REWIRED_PROGRAM" in doc),
    ("parâmetros registrados", "bootstrap_seed" in doc),
]

ok = True
for name, good in checks:
    ok &= good
    print(f"{'PASS' if good else 'FAIL'}  {name}")

sys.exit(0 if ok else 1)
