#!/usr/bin/env bash
# Instala as dependências R do enriquecimento humano (clusterProfiler + org.Hs.eg.db) via BiocManager.
set -euo pipefail

echo "Sistema: $(. /etc/os-release; echo "$PRETTY_NAME")"
echo "R: $(R --version | head -n 1)"

sudo apt-get update
sudo apt-get install -y --no-install-recommends \
    build-essential gfortran cmake \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev \
    libpng-dev libtiff-dev libjpeg-dev libcairo2-dev \
    libglpk-dev libgmp-dev libuv1-dev zlib1g-dev libbz2-dev liblzma-dev libicu-dev

sudo Rscript -e '
options(Ncpus = 2)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
cat("Bioconductor", as.character(BiocManager::version()), "\n")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"), ask = FALSE, update = FALSE)
library(clusterProfiler)
library(org.Hs.eg.db)
cat("OK clusterProfiler", as.character(packageVersion("clusterProfiler")), "\n")
'
