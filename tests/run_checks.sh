#!/usr/bin/env bash
# Roda todos os testes do dataset sintético.
# Uso: bash tests/run_checks.sh [pasta_de_resultados]   (padrão: results/test)
set -euo pipefail

OUT="${1:-results/test}"

python tests/check_truth.py --edges "$OUT/differential_network/differential_edges.tsv" --truth data/example/truth.tsv
python tests/check_topology.py --turnover "$OUT/rewiring/neighbor_turnover.tsv" --metrics "$OUT/topology/network_metrics.tsv"
python tests/check_bootstrap.py --stable "$OUT/rewiring/stable_rewiring.tsv" --edge_stability "$OUT/bootstrap/edge_stability.tsv"
python tests/check_de.py --table "$OUT/rewiring/de_vs_rewiring.tsv"
python tests/check_communities.py --modules "$OUT/communities/gene_modules.tsv"
python tests/check_enrichment.py --enrichment "$OUT/enrichment/custom_enrichment.tsv"
python tests/check_report.py --report "$OUT/report.html"

python tests/check_network_plot.py --png "$OUT/figures/rewiring_network.png" --svg "$OUT/figures/rewiring_network.svg"

echo "Todos os testes passaram."
