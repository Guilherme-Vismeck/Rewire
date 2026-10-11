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


def layout_components(G, seed):
    """Layout de cada componente conexo, empacotados em prateleiras (sem sobreposição)."""
    comps = sorted(nx.connected_components(G), key=lambda c: (-len(c), sorted(c)[0]))
    blocks = []
    for comp in comps:
        n = len(comp)
        nodes = sorted(comp)
        if n == 1:
            local, r = {nodes[0]: np.zeros(2)}, 0.7
        elif n == 2:
            local, r = {nodes[0]: np.array([-0.6, 0.0]), nodes[1]: np.array([0.6, 0.0])}, 1.0
        else:
            r = 1.0 + 0.9 * np.sqrt(n)
            local = nx.spring_layout(G.subgraph(comp), seed=seed, weight="w",
                                     k=1.5 / np.sqrt(n), scale=r, iterations=400)
        blocks.append((local, r))

    gap = 1.2
    total = sum((2 * r + gap) ** 2 for _, r in blocks)
    row_width = max(2 * blocks[0][1] + gap, np.sqrt(total * 1.6))
    pos, x, y, row_h = {}, 0.0, 0.0, 0.0
    for local, r in blocks:
        w = 2 * r + gap
        if x > 0 and x + w > row_width:
            x, y, row_h = 0.0, y - row_h, 0.0
        centre = np.array([x + r, y - r])
        for node, xy in local.items():
            pos[node] = np.asarray(xy) + centre
        x += w
        row_h = max(row_h, w)
    return pos


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
    pos = layout_components(G, args.seed)

    cat = {}
    if de is not None and {"gene", "category"} <= set(de.columns):
        cat = dict(zip(de["gene"].astype(str), de["category"]))
    sizes = {g: 90 + 35 * min(int(counts[g]), 10) for g in genes}
    node_colors = [CATEGORY_COLORS.get(cat.get(g), "#9fb3c8") for g in genes]
    node_sizes = [sizes[g] for g in genes]

    comp_size = {n: len(c) for c in nx.connected_components(G) for n in c}
    if len(genes) <= 15:
        label_genes = list(genes)
    else:
        label_genes = [g for g in genes if counts[g] >= 2 or comp_size[g] <= 3]
    label_genes = label_genes[:args.max_labels]

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
                    if status == "PRESERVED":
                        alpha = 0.55
                    else:
                        alpha = 0.95 if is_stable else 0.45
                    nx.draw_networkx_edges(
                        G, pos, edgelist=list(zip(pp["gene1"], pp["gene2"])),
                        width=(0.6 + 2.6 * pp[rcol].abs()).tolist(),
                        edge_color=STATUS_COLORS[status],
                        style="dashed" if negative else "solid",
                        alpha=alpha, ax=ax)
        nx.draw_networkx_nodes(G, pos, nodelist=genes, node_color=node_colors,
                               node_size=node_sizes, edgecolors="#33415c",
                               linewidths=0.7, ax=ax)
        for g in label_genes:
            r_pt = float(np.sqrt(sizes[g] / np.pi))
            ax.annotate(g, pos[g], xytext=(0, r_pt + 1.5), textcoords="offset points",
                        ha="center", va="bottom", fontsize=7.5, zorder=6,
                        bbox=dict(boxstyle="round,pad=0.12", facecolor="white", edgecolor="none", alpha=0.8))
        ax.set_title(f"{label}  ({len(shown)} edges shown)", fontsize=12)
        ax.set_aspect("equal", adjustable="datalim")
        ax.margins(0.08)
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
    fig.tight_layout(rect=(0, 0.14, 1, 0.95))
    fig.savefig(png, dpi=160, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)
    print(f"Figura: {len(genes)} genes, {len(sub)} arestas candidatas entre eles -> {png}, {svg}")


if __name__ == "__main__":
    main()
