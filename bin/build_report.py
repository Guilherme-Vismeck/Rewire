#!/usr/bin/env python3
"""Relatório HTML autocontido do RewireNF (REPORT)."""
import argparse
import base64
import html
import io
import platform
import subprocess
from datetime import datetime
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

STATUS_ORDER = ["GAINED", "LOST", "SIGN_FLIPPED", "STRENGTHENED", "WEAKENED", "PRESERVED"]
STATUS_COLORS = {"GAINED": "#2a9d8f", "LOST": "#e76f51", "SIGN_FLIPPED": "#9b5de5",
                 "STRENGTHENED": "#457b9d", "WEAKENED": "#e9c46a", "PRESERVED": "#adb5bd"}
CATEGORY_ORDER = ["RW_only", "DEG+RW", "DEG_only", "unchanged"]
CATEGORY_COLORS = {"RW_only": "#e76f51", "DEG+RW": "#9b5de5",
                   "DEG_only": "#457b9d", "unchanged": "#adb5bd"}
LIST_ORDER = ["RW_only", "rewired_stable", "gained_connections", "lost_connections", "module_changed"]

CSS = """
body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;margin:0;color:#1f2933;background:#f7f8fa;line-height:1.5}
main{max-width:980px;margin:0 auto;padding:24px 20px 60px}
h1{margin:0 0 4px;font-size:1.8rem}
h2{margin-top:2.2rem;border-bottom:2px solid #e3e7ec;padding-bottom:6px;font-size:1.25rem}
.sub{color:#52606d;margin:0 0 20px}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px;margin:16px 0}
.card{background:#fff;border:1px solid #e3e7ec;border-radius:8px;padding:12px 14px}
.card b{display:block;font-size:1.5rem}
.card span{color:#52606d;font-size:.85rem}
figure{margin:16px 0;background:#fff;border:1px solid #e3e7ec;border-radius:8px;padding:12px;text-align:center}
figure img{max-width:100%;height:auto}
figcaption{color:#52606d;font-size:.85rem;margin-top:6px}
.tbl{border-collapse:collapse;width:100%;background:#fff;font-size:.85rem}
.tbl th,.tbl td{border-bottom:1px solid #e3e7ec;padding:5px 8px;text-align:left}
.tbl th{background:#eef1f5}
.wrap{overflow-x:auto}
.muted{color:#7b8794;font-size:.85rem}
"""

esc = html.escape


def to_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def read(path, **kwargs):
    p = Path(path)
    if not p.is_file() or p.stat().st_size == 0:
        return None
    try:
        return pd.read_csv(p, sep="\t", **kwargs)
    except Exception:
        return None


def fig_html(fig, caption):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    data = base64.b64encode(buf.getvalue()).decode()
    return (f'<figure><img src="data:image/png;base64,{data}" alt="{esc(caption)}">'
            f"<figcaption>{esc(caption)}</figcaption></figure>")


def table_html(df, max_rows=20):
    if df is None or len(df) == 0:
        return '<p class="muted">No data.</p>'
    d = df.head(max_rows).copy()
    for c in d.select_dtypes(include="float").columns:
        d[c] = d[c].map(lambda v: "" if pd.isna(v) else f"{v:.3g}")
    note = (f'<p class="muted">Showing {len(d)} of {len(df)} rows.</p>' if len(df) > max_rows else "")
    return f'<div class="wrap">{d.to_html(index=False, border=0, classes="tbl", escape=True)}</div>{note}'


def card(label, value):
    return f'<div class="card"><b>{esc(str(value))}</b><span>{esc(label)}</span></div>'


# ---------------- figures ----------------
def fig_correlations(inp, la, lb, min_cor):
    series = {}
    for label in (la, lb):
        p = inp / f"{label}_correlations.tsv.gz"
        if p.is_file():
            r = pd.read_csv(p, sep="\t", usecols=["r"])["r"].to_numpy()
            series[label] = r[~np.isnan(r)]
    if not series:
        return ""
    fig, ax = plt.subplots(figsize=(6.8, 3.6))
    bins = np.linspace(-1, 1, 81)
    for label, r in series.items():
        ax.hist(r, bins=bins, histtype="step", linewidth=1.6, label=f"{label} ({len(r):,} pairs)")
    if min_cor:
        for s in (-min_cor, min_cor):
            ax.axvline(s, color="#666", linestyle=":", linewidth=1)
    ax.set_yscale("log")
    ax.set_xlabel("Correlation (r)")
    ax.set_ylabel("Gene pairs (log scale)")
    ax.legend(frameon=False)
    return fig_html(fig, "Distribution of pairwise correlations in each group (dotted lines: min_cor).")


def fig_status(edges):
    if edges is None or edges.empty:
        return ""
    counts = edges["status"].value_counts().reindex(STATUS_ORDER, fill_value=0)
    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    ax.bar(counts.index, counts.values, color=[STATUS_COLORS[s] for s in counts.index])
    for i, v in enumerate(counts.values):
        ax.text(i, v, str(int(v)), ha="center", va="bottom", fontsize=9)
    ax.set_ylabel("Edges")
    ax.tick_params(axis="x", rotation=25)
    return fig_html(fig, "Candidate edges by status (edges present in at least one network).")


def fig_volcano(edges, fdr, la, lb):
    if edges is None or edges.empty or not {"delta_r", "fdr_diff", "status"} <= set(edges.columns):
        return ""
    y = -np.log10(np.clip(edges["fdr_diff"].astype(float), 1e-50, 1))
    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    for status in STATUS_ORDER:
        m = edges["status"] == status
        if m.any():
            ax.scatter(edges.loc[m, "delta_r"], y[m], s=24, alpha=0.85, color=STATUS_COLORS[status],
                       label=status, edgecolors="none")
    if fdr:
        ax.axhline(-np.log10(fdr), color="#666", linestyle=":", linewidth=1)
    ax.set_xlabel(f"Change in correlation (r in {lb} minus r in {la})")
    ax.set_ylabel("-log10(FDR of the difference)")
    ax.legend(frameon=False, fontsize=8)
    return fig_html(fig, "Differential co-expression of candidate edges (y capped at 50; dotted line: FDR threshold).")


def fig_degree(metrics, la, lb):
    if metrics is None or "delta_degree" not in metrics.columns:
        return ""
    d = metrics["delta_degree"].dropna()
    d = d[d != 0]
    if d.empty:
        return ""
    top = d.reindex(d.abs().sort_values(ascending=False).index).head(15)[::-1]
    fig, ax = plt.subplots(figsize=(6.4, 0.32 * len(top) + 1.3))
    ax.barh(top.index.astype(str), top.values, color=["#2a9d8f" if v > 0 else "#e76f51" for v in top.values])
    ax.set_xlabel(f"Change in degree ({lb} minus {la})")
    return fig_html(fig, "Genes with the largest change in degree between groups.")


def fig_turnover(turn):
    if turn is None or "neighbor_turnover" not in turn.columns:
        return ""
    v = turn["neighbor_turnover"].dropna()
    if v.empty:
        return ""
    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    ax.hist(v, bins=np.linspace(0, 1, 21), color="#457b9d")
    ax.set_xlabel("Neighbor turnover (1 - Jaccard)")
    ax.set_ylabel("Genes")
    return fig_html(fig, f"Neighbor turnover of the {len(v):,} genes with neighbors in at least one network.")


def fig_bootstrap(estab, threshold):
    if estab is None or "bootstrap_support" not in estab.columns or estab.empty:
        return ""
    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    ax.hist(estab["bootstrap_support"].dropna(), bins=np.linspace(0, 1, 21), color="#2a9d8f")
    if threshold:
        ax.axvline(threshold, color="#e76f51", linestyle="--", linewidth=1.2, label=f"threshold {threshold}")
        ax.legend(frameon=False)
    ax.set_xlabel("Bootstrap support")
    ax.set_ylabel("Candidate edges")
    return fig_html(fig, "Bootstrap support of candidate edges.")


def fig_de_rw(de_rw):
    need = {"logFC", "stable_rewired_edges", "category"}
    if de_rw is None or de_rw.empty or not need <= set(de_rw.columns):
        return ""
    rng = np.random.default_rng(1)
    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    for cat in CATEGORY_ORDER:
        m = de_rw["category"] == cat
        if m.any():
            ax.scatter(de_rw.loc[m, "logFC"].abs(),
                       de_rw.loc[m, "stable_rewired_edges"] + rng.uniform(-0.15, 0.15, int(m.sum())),
                       s=22, alpha=0.8, color=CATEGORY_COLORS[cat], label=f"{cat} ({int(m.sum())})",
                       edgecolors="none")
    ax.set_xlabel("|logFC| (differential expression)")
    ax.set_ylabel("Stable rewired edges (jittered)")
    ax.legend(frameon=False, fontsize=8)
    return fig_html(fig, "Differential expression versus stable rewiring, per gene.")


def enrichment_frame(custom, go, kegg):
    frames = []
    if custom is not None and len(custom):
        frames.append(custom.rename(columns={"term": "label"})[["gene_list", "label", "fdr"]].assign(source="custom"))
    for name, df in (("GO", go), ("KEGG", kegg)):
        if df is not None and len(df):
            frames.append(df.rename(columns={"description": "label"})[["gene_list", "label", "fdr"]].assign(source=name))
    return pd.concat(frames) if frames else None


def fig_enrichment(allv):
    if allv is None:
        return ""
    chosen = None
    for gl in LIST_ORDER:
        sub = allv[allv["gene_list"] == gl].sort_values("fdr").head(10)
        if len(sub):
            chosen = (gl, sub)
            break
    if chosen is None:
        return ""
    gl, sub = chosen
    sub = sub[::-1]
    labels = [f"{s[:55]} [{src}]" for s, src in zip(sub["label"].astype(str), sub["source"])]
    fig, ax = plt.subplots(figsize=(7.2, 0.34 * len(sub) + 1.3))
    ax.barh(labels, -np.log10(np.clip(sub["fdr"].astype(float), 1e-300, 1)), color="#9b5de5")
    ax.set_xlabel("-log10(FDR)")
    return fig_html(fig, f"Top enriched terms for the '{gl}' gene list.")


# ---------------- software ----------------
def software_versions():
    rows = [("python", platform.python_version())]
    for pkg in ["pandas", "numpy", "scipy", "statsmodels", "networkx", "matplotlib"]:
        try:
            rows.append((pkg, version(pkg)))
        except PackageNotFoundError:
            pass
    try:
        out = subprocess.run(["R", "--version"], capture_output=True, text=True, timeout=30).stdout
        rows.append(("R", out.splitlines()[0]))
        code = ('for (p in c("limma","clusterProfiler","org.Hs.eg.db")) cat(p, '
                'if (requireNamespace(p, quietly=TRUE)) as.character(packageVersion(p)) '
                'else "not installed", "\\n")')
        out = subprocess.run(["Rscript", "-e", code], capture_output=True, text=True, timeout=60).stdout
        for line in out.splitlines():
            if line.strip():
                name, _, ver = line.partition(" ")
                rows.append((name, ver.strip()))
    except Exception:
        pass
    return pd.DataFrame(rows, columns=["software", "version"])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inputs", default=".")
    ap.add_argument("--parameters", default=None)
    ap.add_argument("--label_a", required=True)
    ap.add_argument("--label_b", required=True)
    ap.add_argument("--pipeline_version", default="unknown")
    ap.add_argument("--nextflow_version", default="unknown")
    ap.add_argument("--output", default="report.html")
    args = ap.parse_args()

    inp = Path(args.inputs)
    la, lb = args.label_a, args.label_b

    params = {}
    if args.parameters and Path(args.parameters).is_file():
        for line in Path(args.parameters).read_text().splitlines():
            if "\t" in line:
                k, v = line.split("\t", 1)
                params[k] = v
    min_cor = to_float(params.get("min_cor"))
    fdr = to_float(params.get("fdr"))
    threshold = to_float(params.get("stability_threshold"))

    qc = read(inp / "qc_report.tsv")
    edges = read(inp / "differential_edges.tsv")
    genes_rw = read(inp / "gene_rewiring.tsv")
    metrics = read(inp / "network_metrics.tsv", index_col=0)
    turn = read(inp / "neighbor_turnover.tsv", index_col=0)
    estab = read(inp / "edge_stability.tsv")
    stable = read(inp / "stable_rewiring.tsv")
    gstab = read(inp / "gene_stability.tsv")
    de_rw = read(inp / "de_vs_rewiring.tsv")
    modules = read(inp / "gene_modules.tsv")
    msum = read(inp / "module_summary.tsv")
    custom = read(inp / "custom_enrichment.tsv")
    go = read(inp / "GO.tsv")
    kegg = read(inp / "KEGG.tsv")
    nets = {label: read(inp / f"{label}_network.tsv") for label in (la, lb)}

    # ---- summary cards ----
    cards = []
    if metrics is not None:
        cards.append(card("genes analyzed", f"{len(metrics):,}"))
    for label in (la, lb):
        if nets[label] is not None:
            cards.append(card(f"edges in {label}", f"{len(nets[label]):,}"))
    if edges is not None and "dcl" in edges.columns:
        cards.append(card("significant differential edges", int(edges["dcl"].astype(bool).sum())))
    if stable is not None:
        cards.append(card("stable rewired edges", len(stable)))
    if de_rw is not None and "category" in de_rw.columns:
        cards.append(card("RW-only genes", int((de_rw["category"] == "RW_only").sum())))
    n_enr = sum(len(d) for d in (custom, go, kegg) if d is not None)
    if any(d is not None for d in (custom, go, kegg)):
        cards.append(card("enriched terms", n_enr))

    body = []

    def section(title, *parts):
        body.append(f"<h2>{esc(title)}</h2>" + "".join(p for p in parts if p))

    section("Dataset and quality control", table_html(qc, 40))

    n_txt = ""
    if all(nets[label] is not None for label in (la, lb)):
        n_txt = (f'<p>Edges present: {esc(la)} {len(nets[la]):,}, {esc(lb)} {len(nets[lb]):,} '
                 f"(|r| &ge; {min_cor}, FDR &lt; {fdr}).</p>")
    section("Co-expression networks", n_txt, fig_correlations(inp, la, lb, min_cor),
            fig_degree(metrics, la, lb))

    sig = None
    if edges is not None and "dcl" in edges.columns:
        cols = [c for c in ["gene1", "gene2", "status", f"r_{la}", f"r_{lb}", "delta_r", "fdr_diff"]
                if c in edges.columns]
        sig = edges[edges["dcl"].astype(bool)][cols]
    section("Differential co-expression", fig_status(edges), fig_volcano(edges, fdr, la, lb),
            "<h3>Significant differential edges</h3>", table_html(sig))

    net_html = ""
    net_png = inp / "rewiring_network.png"
    if net_png.is_file() and net_png.stat().st_size > 0:
        net_b64 = base64.b64encode(net_png.read_bytes()).decode()
        net_html = (
            '<figure><img src="data:image/png;base64,' + net_b64 + '" alt="Rewiring networks">'
            "<figcaption>Co-expression networks of the most rewired genes in each group. "
            "Edge colors show the change in status; dashed edges are negative correlations."
            "</figcaption></figure>"
        )

    top = None
    if genes_rw is not None and not genes_rw.empty:
        top = genes_rw.copy()
        if gstab is not None and len(gstab):
            top = top.merge(gstab[["gene", "stable_rewired_edges", "mean_support"]], on="gene", how="left")
        if metrics is not None and "delta_degree" in metrics.columns:
            top = top.merge(metrics[["delta_degree"]].reset_index(), on="gene", how="left")
        if turn is not None and "neighbor_turnover" in turn.columns:
            top = top.merge(turn[["neighbor_turnover"]].reset_index(), on="gene", how="left")
        top = top.sort_values("dcl_total", ascending=False)
    section("Rewiring and topology", net_html, fig_turnover(turn), "<h3>Top rewired genes</h3>", table_html(top))

    section("Bootstrap stability", fig_bootstrap(estab, threshold),
            "<h3>Stable rewiring edges</h3>", table_html(stable))

    if de_rw is not None:
        summ = (de_rw["category"].value_counts().reindex(CATEGORY_ORDER, fill_value=0)
                .rename_axis("category").reset_index(name="n_genes"))
        rw_cols = [c for c in ["gene", "logFC", "FDR", "stable_rewired_edges", "mean_support"] if c in de_rw.columns]
        rw_only = de_rw[de_rw["category"] == "RW_only"][rw_cols]
        section("Differential expression versus rewiring", fig_de_rw(de_rw), table_html(summ),
                "<h3>RW-only genes</h3>", table_html(rw_only))

    if modules is not None:
        status_counts = (modules["module_status"].value_counts()
                         .rename_axis("module_status").reset_index(name="n_genes"))
        mcols = [c for c in ["gene", f"module_{la}", f"module_{lb}", "module_jaccard", "module_status"]
                 if c in modules.columns]
        changed = modules[modules["module_status"].isin(["changed", "lost_module", "gained_module"])][mcols]
        sizes = None
        if msum is not None and len(msum):
            sizes = msum[["group", "module", "n_genes"]]
        section("Communities", table_html(status_counts), "<h3>Modules</h3>", table_html(sizes),
                "<h3>Genes that change module</h3>", table_html(changed))

    allv = enrichment_frame(custom, go, kegg)
    if any(d is not None for d in (custom, go, kegg)):
        parts = [fig_enrichment(allv)]
        if custom is not None and len(custom):
            cc = custom.sort_values(["gene_list", "fdr"]).groupby("gene_list").head(5)
            parts += ["<h3>Custom gene sets</h3>",
                      table_html(cc[["gene_list", "term", "overlap", "fold_enrichment", "fdr"]], 40)]
        for name, df in (("GO", go), ("KEGG", kegg)):
            if df is not None and len(df):
                dd = df.sort_values(["gene_list", "fdr"]).groupby("gene_list").head(5)
                parts += [f"<h3>{name}</h3>",
                          table_html(dd[["gene_list", "ontology", "description", "count", "fold_enrichment", "fdr"]], 40)]
        section("Functional enrichment", *parts)

    p_df = pd.DataFrame(sorted(params.items()), columns=["parameter", "value"]) if params else None
    section("Parameters and software", "<h3>Parameters</h3>", table_html(p_df, 100),
            "<h3>Software</h3>", table_html(software_versions(), 40))

    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    doc = (
        '<!doctype html><html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width, initial-scale=1">'
        "<title>RewireNF report</title><style>" + CSS + "</style></head><body><main>"
        "<h1>RewireNF report</h1>"
        f'<p class="sub">{esc(la)} vs {esc(lb)} &middot; RewireNF {esc(args.pipeline_version)} '
        f"&middot; Nextflow {esc(args.nextflow_version)} &middot; generated {now}</p>"
        f'<div class="cards">{"".join(cards)}</div>' + "".join(body) + "</main></body></html>"
    )
    Path(args.output).write_text(doc, encoding="utf-8")
    print(f"Relatório escrito em {args.output} ({len(doc) / 1024:.0f} kB)")


if __name__ == "__main__":
    main()
