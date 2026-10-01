#!/usr/bin/env python3
"""Gera um dataset sintético com rewiring conhecido para testar o RewireNF."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def correlated(z, r, rng):
    """Vetor com correlação esperada r com z (z ~ N(0,1))."""
    return r * z + np.sqrt(1 - r**2) * rng.normal(size=z.shape)


def make_group(n, condition, n_background, bg_means, rng):
    is_ctrl = condition == "control"
    genes = {}

    # GENE1 perde GENE2 e ganha GENE3; média igual nas duas condições
    z1 = rng.normal(size=n)
    genes["GENE1"] = 8 + z1
    genes["GENE2"] = 8 + (correlated(z1, 0.9, rng) if is_ctrl else rng.normal(size=n))
    genes["GENE3"] = 8 + (rng.normal(size=n) if is_ctrl else correlated(z1, 0.9, rng))

    # GENE4-GENE5: inversão de sinal
    z4 = rng.normal(size=n)
    genes["GENE4"] = 8 + z4
    genes["GENE5"] = 8 + correlated(z4, 0.9 if is_ctrl else -0.9, rng)

    # GENE6-GENE7: preservado
    z6 = rng.normal(size=n)
    genes["GENE6"] = 8 + z6
    genes["GENE7"] = 8 + correlated(z6, 0.9, rng)

    # GENE8: só expressão diferencial, sem mudar de parceiros
    genes["GENE8"] = (6 if is_ctrl else 8) + rng.normal(size=n)

    # Genes de fundo, independentes
    for i in range(n_background):
        genes[f"BG{i + 1:03d}"] = bg_means[i] + rng.normal(size=n)

    return pd.DataFrame(genes)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--outdir", default="data/example")
    p.add_argument("--n_per_group", type=int, default=40)
    p.add_argument("--n_background", type=int, default=100)
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()

    rng = np.random.default_rng(args.seed)
    bg_means = np.random.default_rng(args.seed + 1).uniform(4, 10, args.n_background)

    ctrl = make_group(args.n_per_group, "control", args.n_background, bg_means, rng)
    case = make_group(args.n_per_group, "case", args.n_background, bg_means, rng)

    ctrl.index = [f"CTRL{i + 1:02d}" for i in range(args.n_per_group)]
    case.index = [f"CASE{i + 1:02d}" for i in range(args.n_per_group)]

    samples_by_genes = pd.concat([ctrl, case])
    expression = samples_by_genes.T.round(3)  # genes nas linhas, amostras nas colunas
    expression.index.name = "gene_id"

    metadata = pd.DataFrame(
        {
            "sample": samples_by_genes.index,
            "condition": ["control"] * args.n_per_group + ["case"] * args.n_per_group,
        }
    )

    truth = pd.DataFrame(
        [
            ("GENE1", "GENE2", "LOST"),
            ("GENE1", "GENE3", "GAINED"),
            ("GENE4", "GENE5", "SIGN_FLIPPED"),
            ("GENE6", "GENE7", "PRESERVED"),
        ],
        columns=["gene_a", "gene_b", "expected"],
    )

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    expression.to_csv(out / "expression.tsv", sep="\t")
    metadata.to_csv(out / "metadata.tsv", sep="\t", index=False)
    truth.to_csv(out / "truth.tsv", sep="\t", index=False)
    print(f"Escrito em {out}: {expression.shape[0]} genes x {expression.shape[1]} amostras")


if __name__ == "__main__":
    main()