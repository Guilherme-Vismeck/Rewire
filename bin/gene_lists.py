#!/usr/bin/env python3
"""Define as listas de genes a enriquecer e o universo (PREPARE_GENE_LISTS)."""
import argparse
from pathlib import Path

import pandas as pd


def endpoints(df):
    return set(df["gene1"].astype(str)) | set(df["gene2"].astype(str))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--de_vs_rewiring", required=True)
    p.add_argument("--stable", required=True, help="stable_rewiring.tsv")
    p.add_argument("--gene_modules", required=True)
    p.add_argument("--outdir", default=".")
    args = p.parse_args()

    table = pd.read_csv(args.de_vs_rewiring, sep="\t")
    stable = pd.read_csv(args.stable, sep="\t")
    modules = pd.read_csv(args.gene_modules, sep="\t")

    universe = table["gene"].astype(str).tolist()
    universe_set = set(universe)

    sets = {
        "rewired_stable": set(table.loc[table["rewired"].astype(bool), "gene"].astype(str)),
        "RW_only": set(table.loc[table["category"] == "RW_only", "gene"].astype(str)),
        "gained_connections": endpoints(stable[stable["status"] == "GAINED"]),
        "lost_connections": endpoints(stable[stable["status"] == "LOST"]),
        "module_changed": set(modules.loc[
            modules["module_status"].isin(["changed", "lost_module", "gained_module"]),
            "gene"].astype(str)),
    }

    rows, sizes = [], []
    for name, genes in sets.items():
        genes = sorted(genes & universe_set)
        sizes.append((name, len(genes)))
        rows += [(name, g) for g in genes]

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=["set", "gene"]).to_csv(out / "gene_lists.tsv", sep="\t", index=False)
    pd.DataFrame(sizes, columns=["gene_list", "n_genes"]).to_csv(
        out / "gene_list_sizes.tsv", sep="\t", index=False)
    (out / "universe.txt").write_text("\n".join(universe) + "\n")

    print(f"Universo: {len(universe)} genes")
    for name, n in sizes:
        print(f"  {name}: {n} genes")


if __name__ == "__main__":
    main()
