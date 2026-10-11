#!/usr/bin/env python3
"""Confere a figura das redes no dataset sintético."""
import argparse
import sys
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument("--png", required=True)
p.add_argument("--svg", required=True)
args = p.parse_args()

png = Path(args.png).read_bytes()
svg = Path(args.svg).read_text(encoding="utf-8")

checks = [
    ("PNG válido", png[:8] == b"\x89PNG\r\n\x1a\n"),
    ("PNG com tamanho razoável (> 20 kB)", len(png) > 20_000),
    ("SVG contém os genes rewired", all(g in svg for g in ["GENE1", "GENE2", "GENE3", "GENE4", "GENE5"])),
    ("SVG não contém genes sem rewiring", not any(g in svg for g in ["GENE6", "GENE7", "GENE8", "BG001"])),
]

ok = True
for name, good in checks:
    ok &= good
    print(f"{'PASS' if good else 'FAIL'}  {name}")

sys.exit(0 if ok else 1)
