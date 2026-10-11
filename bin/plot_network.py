#!/usr/bin/env python3
"""Figura das redes de coexpressão dos genes mais rewired, uma por grupo (PLOT_NETWORK)."""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import networkx as nx  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

plt.rcParams["svg.fonttype"] = "none"  # mantém os nomes dos genes como texto no SVG

STATUS_COLORS = {"GAINED": "#2a9d8f", "LOST": "#e76f51", "SIGN_FLIPPED": "#9b5de5",
                 "STRENGTHENED": "#457b9d", "WEAKENED": "#e9c46a", "PRESERVED": "#b8c0c8"}
DRAW_ORDER = ["PRESERVED", "WEAKENED", "STRENGTHENED", "SIGN_FLIPPED", "GAINED", "LOST"]
CHANGED = ["GAINED", "LOST", "SIGN_FLIPPED", "STRENGTHENED", "WEAKENED"]
CATEGORY_COLORS = {"RW_only": "#e76f51", "DEG+RW": "#9b5de5",
                   "DEG_only": "#457b9d", "unchanged": "#dfe3e8"}


def read_optional(path):
    p = Path(path) if path else None
    if p is None or not p.is_file() or p.stat().st_size == 0:
        return None
    try:
        return pd.read_csv(p, sep="\t")
    except Exception:
        return None


def message_figure(png, svg, text):
    fig, ax = plt.subplots(figsize=(7, 2.5))
    ax.text(0.5, 0.5, text, ha="center", va="center", fontsize=12)
    ax.set_axis_off()
    fig.savefig(png, dpi=130, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)
    print(text)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--edges", required=True, help="differential_edges.tsv")
    ap.add_argument("--stable", default=None, help="stable_rewiring.tsv")
    ap.add_argument("--de_vs_rewiring", default=None, help="de_vs_rewiring.tsv (opcional)")
    ap.add_argument("--label_a", required=True)
    ap.add_argument("--label_b", required=True)
    ap.add_argument("--max_genes", type=int, default=40)
    ap.add_argument("--max_labels", type=int, default=30)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--output_prefix", default="rewiring_network")
    args = ap.parse_args()

    la, lb = args.label_a, args.label_b
    png, svg = f"{args.output_prefix}.png", f"{args.output_prefix}.svg"

    edges = read_optional(args.edges)
    stable = read_optional(args.stable)
    de = read_optional(args.de_vs_rewiring)

    if edges is None or edges.empty:
        message_figure(png, svg, "No candidate edges to display")
        return

    edges = edges.copy()
    edges["gene1"] = edges["gene1"].astype(str)
    edges["gene2"] = edges["gene2"].astype(str)
    edges["dcl"] = edges["dcl"].fillna(False).astype(bool)
    pa, pb, ra, rb = f"present_{la}", f"present_{lb}", f"r_{la}", f"r_{lb}"
    for col in (pa, pb):
        edges[col] = edges[col].fillna(False).astype(bool)

    stable_pairs = set()
    if stable is not None and len(stable):
        source = stable
        stable_pairs = {frozenset(p) for p in zip(stable["gene1"].astype(str),
                                                  stable["gene2"].astype(str))}
        basis = "stable rewiring"
    else:
        source = edges[edges["dcl"] & edges["status"].isin(CHANGED)]
        basis = "raw rewiring, no stable edges"
    if source.empty:
        message_figure(png, svg, "No rewired edges to display")
        return

    counts = pd.concat([source["gene1"], source["gene2"]]).astype(str).value_counts()
    genes = list(counts.head(args.max_genes).index)
    gene_set = set(genes)

    sub = edges[edges["gene1"].isin(gene_set) & edges["gene2"].isin(gene_set)].copy()
    if stable_pairs:
        sub["is_stable"] = [frozenset(p) in stable_pairs for p in zip(sub["gene1"], sub["gene2"])]
    else:
        sub["is_stable"] = True

    G = nx.Graph()
    G.add_nodes_from(genes)
    for g1, g2, wa, wb in zip(sub["gene1"], sub["gene2"], sub[ra].abs(), sub[rb].abs()):
        G.add_edge(g1, g2, w=max(float(wa), float(wb)))
    pos = nx.spring_layout(G, seed=args.seed, weight="w",
                           k=1.6 / np.sqrt(max(len(G), 1)), iterations=400)

    cat = {}
    if de is not None and {"gene", "category"} <= set(de.columns):
        cat = dict(zip(de["gene"].astype(str), de["category"]))
    node_colors = [CATEGORY_COLORS.get(cat.get(g), "#9fb3c8") for g in genes]
    node_sizes = [140 + 70 * min(int(counts[g]), 12) for g in genes]
    labels = {g: g for g in genes[:args.max_labels]}

    fig, axes = plt.subplots(1, 2, figsize=(14, 7.4))
    for ax, label, rcol, pcol in ((axes[0], la, ra, pa), (axes[1], lb, rb, pb)):
        shown = sub[sub[pcol]]
        for status in DRAW_ORDER:
            for negative in (False, True):
                part = shown[(shown["status"] == status) & ((shown[rcol] < 0) == negative)]
                for is_stable in (False, True):
                    pp = part[part["is_stable"] == is_stable]
                    if pp.empty:
                        continue
                    faint = status in CHANGED and not is_stable
                    nx.draw_networkx_edges(
                        G, pos, edgelist=list(zip(pp["gene1"], pp["gene2"])),
                        width=(0.8 + 3.2 * pp[rcol].abs()).tolist(),
                        edge_color=STATUS_COLORS[status],
                        style="dashed" if negative else "solid",
                        alpha=0.4 if faint else 0.9, ax=ax)
        nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes,
                               edgecolors="#33415c", linewidths=0.7, ax=ax)
        nx.draw_networkx_labels(G, pos, labels=labels, font_size=8, ax=ax)
        ax.set_title(f"{label}  ({len(shown)} edges shown)", fontsize=12)
        ax.margins(0.12)
        ax.set_axis_off()

    handles = [Line2D([0], [0], color=STATUS_COLORS[s], lw=3, label=t) for s, t in [
        ("LOST", f"lost (only in {la})"), ("GAINED", f"gained (only in {lb})"),
        ("SIGN_FLIPPED", "sign flipped"), ("STRENGTHENED", "strengthened"),
        ("WEAKENED", "weakened"), ("PRESERVED", "preserved")]]
    handles += [Line2D([0], [0], color="#555", lw=2, linestyle="solid", label="positive correlation"),
                Line2D([0], [0], color="#555", lw=2, linestyle="dashed", label="negative correlation")]
    if cat:
        handles += [Line2D([0], [0], marker="o", linestyle="", markerfacecolor=c,
                           markeredgecolor="#33415c", markersize=9, label=name)
                    for name, c in CATEGORY_COLORS.items()]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=9)
    fig.suptitle(f"Co-expression networks of the {len(genes)} most rewired genes ({basis}); "
                 "node size = number of rewired edges", fontsize=11)
    fig.tight_layout(rect=(0, 0.16, 1, 0.95))
    fig.savefig(png, dpi=160, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)
    print(f"Figura: {len(genes)} genes, {len(sub)} arestas candidatas entre eles -> {png}, {svg}")


if __name__ == "__main__":
    main()
