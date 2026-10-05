#!/usr/bin/env bash
# Testa enrich_human.R (GO) com genes reais de ciclo celular, sem internet.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
tmp="$(mktemp -d)"
cd "$tmp"

Rscript -e '
suppressPackageStartupMessages(library(org.Hs.eg.db))
cc <- c("CDK1","CCNB1","CCNB2","CDC20","PLK1","AURKA","AURKB","BUB1","BUB1B",
        "MAD2L1","CCNA2","CDC25C","TOP2A","KIF11","CENPE")
set.seed(1)
bg <- sample(keys(org.Hs.eg.db, keytype = "SYMBOL"), 3000)
writeLines(unique(c(cc, bg)), "universe.txt")
write.table(data.frame(set = "cell_cycle", gene = cc), "gene_lists.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
'

"$REPO/bin/enrich_human.R" --gene_lists gene_lists.tsv --universe universe.txt \
    --gene_id_type SYMBOL --run_kegg false --fdr 0.05

python3 - "$tmp/GO.tsv" << 'PY'
import sys
import pandas as pd

go = pd.read_csv(sys.argv[1], sep="\t")
hits = go[(go["ontology"] == "GO_BP") &
          go["description"].str.contains("mitotic|cell cycle|division|segregation",
                                         case=False, regex=True)]
print(f"Termos GO BP de ciclo celular/divisão significativos: {len(hits)}")
print(hits[["description", "count", "fdr"]].head(5).to_string(index=False))
sys.exit(0 if len(hits) > 0 else 1)
PY
echo "PASS  enrich_human.R: GO de ciclo celular encontrado"
