FROM rocker/r-ver:4.5.2

LABEL org.opencontainers.image.title="RewireNF" \
      org.opencontainers.image.description="Runtime environment for the RewireNF pipeline" \
      org.opencontainers.image.source="https://github.com/Guilherme-Vismeck/Rewire" \
      org.opencontainers.image.licenses="MIT"

ENV DEBIAN_FRONTEND=noninteractive \
    HOME=/tmp \
    MPLCONFIGDIR=/tmp/matplotlib

# System libraries: build tools for compiled R packages, plus procps (required by Nextflow)
RUN apt-get update && apt-get install -y --no-install-recommends \
        python3 python3-venv python3-pip procps \
        build-essential gfortran cmake \
        libcurl4-openssl-dev libssl-dev libxml2-dev libuv1-dev \
        libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev \
        libpng-dev libtiff-dev libjpeg-dev libcairo2-dev \
        libglpk-dev libgmp-dev zlib1g-dev libbz2-dev liblzma-dev libicu-dev \
    && rm -rf /var/lib/apt/lists/*

# Python environment (versions pinned in requirements.txt)
COPY requirements.txt /tmp/requirements.txt
RUN python3 -m venv /opt/venv \
    && /opt/venv/bin/pip install --no-cache-dir -r /tmp/requirements.txt \
    && rm /tmp/requirements.txt
ENV PATH="/opt/venv/bin:${PATH}"

# R packages: the build fails if any of them is missing
RUN Rscript -e ' \
      options(Ncpus = 2); \
      install.packages("BiocManager"); \
      pkgs <- c("limma", "clusterProfiler", "org.Hs.eg.db"); \
      BiocManager::install(pkgs, ask = FALSE, update = FALSE); \
      missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]; \
      if (length(missing)) stop("Failed to install: ", paste(missing, collapse = ", ")) \
    '

CMD ["bash"]
