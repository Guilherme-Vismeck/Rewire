#!/usr/bin/env python3
"""Dataset de demonstração: um gene hub troca de módulo entre as condições."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def make_group(n, is_control, n_background, bg_means, rng, rho=0.8):
    def load(f):
        return np.sqrt(rho) * f + np.sqrt(1 - rho) * rng.normal(size=n)

    def factor():
        f = rng.normal(size=n)
        return (f - f.mean()) / f.std()  # sem deslocamento de média entre os grupos

    fa, fb, fc = factor(), factor(), factor()
    genes = {}

    # Três módulos que existem nas duas condições
    for i in range(1, 9):
        genes[f"MODA_{i}"] = 8 + load(fa)
    for i in range(1, 9):
        genes[f"MODB_{i}"] = 8 + load(fb)
    for i in range(1, 7):
        genes[f"MODC_{i}"] = 8 + load(fc)

    # HUB1: módulo A no controle, módulo B no caso (mesma expressão média)
    genes["HUB1"] = 8 + load(fa if is_control else fb)
    # HUB2: módulo C no controle, nenhum parceiro no caso
    genes["HUB2"] = 8 + (load(fc) if is_control else rng.normal(size=n))

    # Inversão de sinal
    z = rng.normal(size=n)
    genes["FLIP_1"] = 8 + z
    genes["FLIP_2"] = 8 + (0.9 if is_control else -0.9) * z + np.sqrt(1 - 0.81) * rng.normal(size=n)

    # Expressão diferencial sem rewiring
    for i in range(1, 7):
        genes[f"DEG_{i}"] = (6 if is_control else 8) + rng.normal(size=n)

    for i in range(n_background):
        genes[f"BG{i + 1:03d}"] = bg_means[i] + rng.normal(size=n)
    return pd.DataFrame(genes)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--outdir", default="data/demo")
    p.add_argument("--n_per_group", type=int, default=60)
    p.add_argument("--n_background", type=int, default=60)
    p.add_argument("--seed", type=int, default=7)
    args = p.parse_args()

    rng = np.random.default_rng(args.seed)
    bg_means = np.random.default_rng(args.seed + 1).uniform(4, 10, args.n_background)

    ctrl = make_group(args.n_per_group, True, args.n_background, bg_means, rng)
    case = make_group(args.n_per_group, False, args.n_background, bg_means, rng)
    ctrl.index = [f"CTRL{i + 1:02d}" for i in range(args.n_per_group)]
    case.index = [f"CASE{i + 1:02d}" for i in range(args.n_per_group)]

    both = pd.concat([ctrl, case])
    expression = both.T.round(3)
    expression.index.name = "gene_id"
    metadata = pd.DataFrame({
        "sample": both.index,
        "condition": ["control"] * args.n_per_group + ["case"] * args.n_per_group,
    })

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    expression.to_csv(out / "expression.tsv", sep="\t")
    metadata.to_csv(out / "metadata.tsv", sep="\t", index=False)
    print(f"Escrito em {out}: {expression.shape[0]} genes x {expression.shape[1]} amostras")


if __name__ == "__main__":
    main()
