# Rewire

![test](https://github.com/Guilherme-Vismeck/Rewire/actions/workflows/test.yml/badge.svg)

A modular Nextflow pipeline for identifying differential gene co-expression and network rewiring between two biological conditions.

RewireNF asks which genes change their *relationships* with other genes between conditions, even when their mean expression stays the same.

> **Status:** v1.0.

## What it does (v0.1)

Starting from a normalized expression matrix and sample metadata:

1. **VALIDATE_INPUTS**: checks files, sample matching, duplicates, missing/infinite values, zero variance and group sizes (warns when n < 30).
2. **PREPROCESS**: removes genes with too many missing values or zero variance, imputes remaining NAs, keeps the most variable genes (variance averaged within groups, so differentially expressed genes are not favored) and splits the groups.
3. **BUILD_NETWORK**: computes correlations, p-values and BH-FDR for every gene pair in each group; an edge is present when |r| >= `min_cor` and FDR < `fdr`.
4. **DIFFERENTIAL_NETWORK**: compares the two networks with the Fisher z-test and classifies each edge as `PRESERVED`, `GAINED`, `LOST`, `STRENGTHENED`, `WEAKENED` or `SIGN_FLIPPED`; summarizes changes per gene.
5. **TOPOLOGY**: computes per-gene degree, strength, betweenness, closeness, eigenvector centrality and clustering in each network (with deltas between groups) and the neighbor turnover (1 - Jaccard similarity of each gene's neighbors).
6. **BOOTSTRAP**: resamples the samples of each group with replacement and measures how often each candidate edge reappears (edge stability). Combined with the differential test, it separates raw rewiring from stable rewiring (`bootstrap_support >= stability_threshold`).
7. **DIFFERENTIAL_EXPRESSION** and **INTEGRATE_DE_RW**: runs limma on the analyzed genes and crosses the result with stable rewiring, classifying each gene as `DEG+RW`, `RW_only` (rewired without differential expression), `DEG_only` or `unchanged`.
8. **COMMUNITIES**: detects modules (Louvain) in each network and flags genes whose module changes between conditions.
9. **ENRICHMENT** (optional): functional enrichment of the rewired gene lists, with your own gene sets and/or human GO/KEGG.
10. **EXPORT_NETWORK**: writes GraphML files for Cytoscape.
11. **REPORT**: builds a self-contained `report.html` (summary cards, figures and tables for every step).

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
| `--community_resolution` | `1.0` | Louvain resolution |
| `--community_seed` | `42` | Random seed for Louvain |
| `--module_overlap_threshold` | `0.5` | Minimum Jaccard overlap to call a gene's module unchanged |
| `--run_enrichment` | `false` | Run functional enrichment |
| `--organism` | none | `human` enables GO/KEGG |
| `--gene_sets` | none | Your own gene sets (TSV with `term`, `gene`, or `.gmt`) |
| `--gene_id_type` | `SYMBOL` | `SYMBOL`, `ENSEMBL` or `ENTREZID` (human) |
| `--run_kegg` | `true` | Include KEGG in the human enrichment |
| `--enrich_fdr` | `0.05` | FDR threshold for enrichment |
| `--enrich_min_size`, `--enrich_max_size` | `5`, `500` | Gene set size limits |
| `--enrich_min_overlap` | `2` | Minimum genes in common to report a term |
| `--run_report` | `true` | Build the HTML report |
| `--container` | `rewirenf:latest` | Image used with `-profile docker` |
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
├── communities/               # gene_modules, module_summary
├── enrichment/                # gene_lists, custom_enrichment, GO, KEGG
├── bootstrap/                 # edge_stability, gene_stability
├── topology/network_metrics.tsv
├── report.html
└── networks_graphml/          # <group_a>.graphml, <group_b>.graphml, rewiring.graphml
```

`rewiring.graphml` carries per-gene `degree_<group>`, `delta_degree`, `edges_gained`, `edges_lost` and `sign_flips`, and per-edge `status`, `delta_r` and `fdr_diff`.

## Functional enrichment

Enabled with `--run_enrichment` (requires `--run_de true`). Five gene lists are tested for over-representation (hypergeometric test, BH-FDR) against a universe made of the genes analyzed by the pipeline: `rewired_stable`, `RW_only`, `gained_connections`, `lost_connections` and `module_changed`.

- **Your own gene sets (any organism):** `--gene_sets sets.tsv`, a TSV with columns `term` and `gene` (or a `.gmt` file).
- **Human GO and KEGG:** `--organism human` runs clusterProfiler (GO BP/MF/CC and KEGG). Set `--gene_id_type` to `SYMBOL` (default), `ENSEMBL` or `ENTREZID`. KEGG queries the KEGG web service; if it is unreachable the run continues with an empty KEGG table (`--run_kegg false` skips it).

```bash
nextflow run main.nf ... --run_enrichment --organism human --gene_sets my_sets.tsv
```

The human path needs R packages that are installed from Bioconductor (compiled, can take a while):

```bash
bash scripts/install_r_deps.sh
```

## Running with Docker

The `Dockerfile` pins Python, R, limma, clusterProfiler and org.Hs.eg.db, so the whole pipeline (including human GO/KEGG enrichment) runs without installing anything else:

```bash
docker build -t rewirenf:latest .
nextflow run main.nf -profile test,docker
```

Use `--container <image>` to run a different image.

## Testing

The synthetic dataset has known ground truth (`data/example/truth.tsv`): GENE1-GENE2 is lost, GENE1-GENE3 is gained, GENE4-GENE5 flips sign and GENE6-GENE7 is preserved. GitHub Actions runs the pipeline on every push and checks these results with `tests/check_truth.py`.

Human GO enrichment is tested locally with `tests/test_enrich_human.sh` (requires the R packages above); the CI runs the custom gene-set enrichment test.

The `docker` workflow builds the image and runs the whole pipeline plus all checks (including human GO enrichment) inside it.

## Development environment

The repository includes a dev container (`.devcontainer/`) with Java, Python, R (limma) and Nextflow, so it works in GitHub Codespaces without local installation. Python dependencies: pandas, numpy, scipy, statsmodels, networkx.

## Limitations

- Starts from a processed expression matrix, not from FASTQ files.
- The Fisher z-test is exact for Pearson correlation and only approximate with Spearman.
- Correlation networks estimated from few samples are unstable; the validation step warns when a group has fewer than 30 samples.
- Bootstrap stability is computed only for candidate edges (present in at least one network) and uses the |r| threshold without FDR, because resampling with replacement makes p-values over-optimistic.
- Built-in functional annotation exists only for human; other organisms need `--gene_sets`.

## Roadmap

- v0.2 (done): topology metrics, neighbor turnover
- v0.3 (done): bootstrap edge stability, stable rewiring
- v0.4 (done): differential expression (limma), DEG vs rewiring categories
- v0.5 (done): communities, GO/KEGG enrichment, custom gene-set enrichment
- v0.6 (done): self-contained HTML report
- v0.7 (done): Docker image and containerized CI
- v1.0 (done): final documentation (citation file, example output)

## Example output

`docs/example/` holds the output of the test run on the synthetic dataset: the HTML report (`report.html`, download it and open it in a browser), the differential edges, the stable rewiring edges, the DE versus rewiring table, the custom enrichment results and the `rewiring.graphml` file for Cytoscape.

## Citation

If you use RewireNF, please cite it with the metadata in `CITATION.cff` (GitHub shows it under "Cite this repository").

## License

MIT
