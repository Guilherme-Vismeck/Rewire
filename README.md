# Rewire

![test](https://github.com/Guilherme-Vismeck/Rewire/actions/workflows/test.yml/badge.svg)

A modular Nextflow pipeline for identifying differential gene co-expression and network rewiring between two biological conditions.

RewireNF asks which genes change their *relationships* with other genes between conditions, even when their mean expression stays the same.

> **Status:** v0.4. Communities and functional enrichment are planned (see Roadmap).

## What it does (v0.1)

Starting from a normalized expression matrix and sample metadata:

1. **VALIDATE_INPUTS**: checks files, sample matching, duplicates, missing/infinite values, zero variance and group sizes (warns when n < 30).
2. **PREPROCESS**: removes genes with too many missing values or zero variance, imputes remaining NAs, keeps the most variable genes (variance averaged within groups, so differentially expressed genes are not favored) and splits the groups.
3. **BUILD_NETWORK**: computes correlations, p-values and BH-FDR for every gene pair in each group; an edge is present when |r| >= `min_cor` and FDR < `fdr`.
4. **DIFFERENTIAL_NETWORK**: compares the two networks with the Fisher z-test and classifies each edge as `PRESERVED`, `GAINED`, `LOST`, `STRENGTHENED`, `WEAKENED` or `SIGN_FLIPPED`; summarizes changes per gene.
5. **TOPOLOGY**: computes per-gene degree, strength, betweenness, closeness, eigenvector centrality and clustering in each network (with deltas between groups) and the neighbor turnover (1 - Jaccard similarity of each gene's neighbors).
6. **BOOTSTRAP**: resamples the samples of each group with replacement and measures how often each candidate edge reappears (edge stability). Combined with the differential test, it separates raw rewiring from stable rewiring (`bootstrap_support >= stability_threshold`).
7. **DIFFERENTIAL_EXPRESSION** and **INTEGRATE_DE_RW**: runs limma on the analyzed genes and crosses the result with stable rewiring, classifying each gene as `DEG+RW`, `RW_only` (rewired without differential expression), `DEG_only` or `unchanged`.
8. **EXPORT_NETWORK**: writes GraphML files for Cytoscape.

## Input

- `expression.tsv`: genes in rows, samples in columns, first column is the gene ID. Expects continuous normalized values (log2(TPM+1), logCPM, VST), not raw counts.
- `metadata.tsv`: one row per sample, with a `sample` column and a grouping column (default `condition`).

## Quick start

```bash
nextflow run main.nf -profile test
```

This runs the pipeline on a small synthetic dataset (`data/example/`) generated with `bin/make_example_data.py`, where the expected result is known.

On your own data:

```bash
nextflow run main.nf \
    --expression data/expression.tsv \
    --metadata data/metadata.tsv \
    --group_col condition \
    --group_a control \
    --group_b disease \
    --method pearson \
    --min_cor 0.6 \
    --fdr 0.05 \
    --outdir results
```

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--expression` | required | Expression matrix (TSV) |
| `--metadata` | required | Sample metadata (TSV) |
| `--group_col` | `condition` | Metadata column with the groups |
| `--group_a`, `--group_b` | required | The two groups to compare |
| `--method` | `pearson` | `pearson` or `spearman` |
| `--min_cor` | `0.6` | Minimum \|r\| for an edge |
| `--fdr` | `0.05` | FDR threshold |
| `--top_variable_genes` | `3000` | Number of most variable genes kept |
| `--bootstrap` | `100` | Number of bootstrap replicates |
| `--bootstrap_seed` | `42` | Random seed for the bootstrap |
| `--stability_threshold` | `0.8` | Minimum bootstrap support for stable rewiring |
| `--run_de` | `true` | Run limma and the DE vs rewiring integration |
| `--de_fdr` | `0.05` | FDR threshold for differential expression |
| `--min_logfc` | `0` | Minimum \|logFC\| to call a DEG |
| `--min_rewired_edges` | `1` | Stable rewired edges needed to call a gene rewired |
| `--outdir` | `results` | Output directory |

## Output

```
results/
├── qc/qc_report.tsv
├── preprocessing/filtered_expression.tsv
├── networks/                  # <group>_network.tsv
├── differential_network/      # differential_edges, gained_edges, lost_edges, sign_flips, gene_rewiring
├── rewiring/neighbor_turnover.tsv
├── rewiring/stable_rewiring.tsv
├── rewiring/de_vs_rewiring.tsv
├── differential_expression/DE_results.tsv
├── bootstrap/                 # edge_stability, gene_stability
├── topology/network_metrics.tsv
└── networks_graphml/          # <group_a>.graphml, <group_b>.graphml, rewiring.graphml
```

`rewiring.graphml` carries per-gene `degree_<group>`, `delta_degree`, `edges_gained`, `edges_lost` and `sign_flips`, and per-edge `status`, `delta_r` and `fdr_diff`.

## Testing

The synthetic dataset has known ground truth (`data/example/truth.tsv`): GENE1-GENE2 is lost, GENE1-GENE3 is gained, GENE4-GENE5 flips sign and GENE6-GENE7 is preserved. GitHub Actions runs the pipeline on every push and checks these results with `tests/check_truth.py`.

## Development environment

The repository includes a dev container (`.devcontainer/`) with Java, Python, R (limma) and Nextflow, so it works in GitHub Codespaces without local installation. Python dependencies: pandas, numpy, scipy, statsmodels, networkx.

## Limitations

- Starts from a processed expression matrix, not from FASTQ files.
- The Fisher z-test is exact for Pearson correlation and only approximate with Spearman.
- Correlation networks estimated from few samples are unstable; the validation step warns when a group has fewer than 30 samples.
- Bootstrap stability is computed only for candidate edges (present in at least one network) and uses the |r| threshold without FDR, because resampling with replacement makes p-values over-optimistic.

## Roadmap

- v0.2 (done): topology metrics, neighbor turnover
- v0.3 (done): bootstrap edge stability, stable rewiring
- v0.4 (done): differential expression (limma), DEG vs rewiring categories
- v0.5: communities, GO/KEGG enrichment
- v1.0: HTML report, Docker, documentation

## License

MIT
