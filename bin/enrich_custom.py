#!/usr/bin/env python3
"""Enriquecimento hipergeométrico com conjuntos de genes do usuário (ENRICH_CUSTOM)."""
import argparse
from pathlib import Path

import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

COLUMNS = ["gene_list", "term", "set_size", "list_size", "overlap",
           "fold_enrichment", "pvalue", "fdr", "genes"]


def read_gene_sets(path):
    """Aceita TSV longo (colunas term e gene) ou GMT."""
    sets = {}
    if path.lower().endswith(".gmt"):
        for line in Path(path).read_text().splitlines():
            parts = line.split("\t")
            if len(parts) >= 3:
                sets.setdefault(parts[0], set()).update(g.strip() for g in parts[2:] if g.strip())
    else:
        df = pd.read_csv(path, sep="\t", dtype=str)
        if not {"term", "gene"} <= set(df.columns):
            raise SystemExit("ERRO: o arquivo de conjuntos precisa das colunas 'term' e 'gene' (ou ser .gmt)")
        for term, gene in zip(df["term"], df["gene"]):
            if isinstance(gene, str) and gene.strip():
                sets.setdefault(term, set()).add(gene.strip())
    return sets


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gene_lists", required=True)
    p.add_argument("--universe", required=True)
    p.add_argument("--gene_sets", required=True)
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--min_size", type=int, default=5)
    p.add_argument("--max_size", type=int, default=500)
    p.add_argument("--min_overlap", type=int, default=2)
    p.add_argument("--outdir", default=".")
    args = p.parse_args()

    universe = {g.strip() for g in Path(args.universe).read_text().splitlines() if g.strip()}
    N = len(universe)
    if N == 0:
        raise SystemExit("ERRO: universo de genes vazio")

    all_sets = read_gene_sets(args.gene_sets)
    sets = {}
    for term, members in all_sets.items():
        inside = members & universe
        if args.min_size <= len(inside) <= args.max_size:
            sets[term] = inside
    print(f"Conjuntos lidos: {len(all_sets)}; usados (tamanho {args.min_size}-{args.max_size} "
          f"dentro do universo): {len(sets)}")
    if not sets:
        print("AVISO: nenhum conjunto de genes com tamanho válido dentro do universo")

    lists = pd.read_csv(args.gene_lists, sep="\t", dtype=str)
    rows = []
    for name, grp in lists.groupby("set", sort=False):
        genes = set(grp["gene"]) & universe
        n = len(genes)
        if n == 0:
            continue
        tested = []
        for term, members in sets.items():
            overlap = genes & members
            k = len(overlap)
            if k == 0:
                continue
            K = len(members)
            tested.append({
                "gene_list": name, "term": term, "set_size": K, "list_size": n,
                "overlap": k, "fold_enrichment": (k / n) / (K / N),
                "pvalue": stats.hypergeom.sf(k - 1, N, K, n),
                "genes": ",".join(sorted(overlap)),
            })
        if not tested:
            continue
        fdr = multipletests([t["pvalue"] for t in tested], method="fdr_bh")[1]
        for t, q in zip(tested, fdr):
            t["fdr"] = q
            if q < args.fdr and t["overlap"] >= args.min_overlap:
                rows.append(t)

    df = pd.DataFrame(rows, columns=COLUMNS)
    if len(df):
        df = df.sort_values(["gene_list", "fdr"])

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "custom_enrichment.tsv", sep="\t", index=False, float_format="%.4g")
    print(f"Termos significativos (FDR < {args.fdr}, sobreposição >= {args.min_overlap}): {len(df)}")


if __name__ == "__main__":
    main()
