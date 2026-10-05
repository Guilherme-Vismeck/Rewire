#!/usr/bin/env Rscript
# Enriquecimento GO/KEGG em humano com clusterProfiler (ENRICH_HUMAN)
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  i <- which(args == paste0("--", name))
  if (length(i) == 0) return(default)
  args[i[1] + 1]
}

gene_lists_file <- get_arg("gene_lists")
universe_file   <- get_arg("universe")
key_type        <- toupper(get_arg("gene_id_type", "SYMBOL"))
run_kegg        <- tolower(get_arg("run_kegg", "true")) == "true"
fdr_cutoff      <- as.numeric(get_arg("fdr", "0.05"))
min_size        <- as.integer(get_arg("min_size", "5"))
max_size        <- as.integer(get_arg("max_size", "500"))
min_overlap     <- as.integer(get_arg("min_overlap", "2"))

if (is.null(gene_lists_file) || is.null(universe_file)) {
  stop("Uso: enrich_human.R --gene_lists X --universe Y [--gene_id_type SYMBOL] [--run_kegg true]")
}
if (!key_type %in% c("SYMBOL", "ENSEMBL", "ENTREZID")) {
  stop("--gene_id_type deve ser SYMBOL, ENSEMBL ou ENTREZID")
}

clean_ids <- function(x) {
  x <- unique(trimws(as.character(x)))
  x <- x[nzchar(x)]
  if (key_type == "ENSEMBL") x <- sub("\\..*$", "", x)  # remove versão (ENSG...13)
  unique(x)
}

gene_lists <- read.delim(gene_lists_file, stringsAsFactors = FALSE)
universe <- clean_ids(readLines(universe_file))
message("Universo: ", length(universe), " genes (", key_type, ")")

empty_result <- data.frame(
  gene_list = character(), ontology = character(), term_id = character(),
  description = character(), count = integer(), fold_enrichment = numeric(),
  pvalue = numeric(), fdr = numeric(), genes = character(), stringsAsFactors = FALSE
)

ratio <- function(x) {
  vapply(strsplit(x, "/"), function(p) as.numeric(p[1]) / as.numeric(p[2]), numeric(1))
}

format_result <- function(res, set, ontology, id_map = NULL) {
  if (is.null(res)) return(NULL)
  df <- as.data.frame(res)
  if (nrow(df) == 0) return(NULL)
  df <- df[df$p.adjust < fdr_cutoff & df$Count >= min_overlap, , drop = FALSE]
  if (nrow(df) == 0) return(NULL)
  genes <- if (is.null(id_map)) {
    gsub("/", ",", df$geneID)
  } else {
    vapply(strsplit(df$geneID, "/"),
           function(g) paste(unique(id_map[g]), collapse = ","), character(1))
  }
  data.frame(
    gene_list = set, ontology = ontology, term_id = df$ID, description = df$Description,
    count = df$Count, fold_enrichment = ratio(df$GeneRatio) / ratio(df$BgRatio),
    pvalue = df$pvalue, fdr = df$p.adjust, genes = genes, stringsAsFactors = FALSE
  )
}

list_genes <- function(set) {
  intersect(clean_ids(gene_lists$gene[gene_lists$set == set]), universe)
}

# ---------------- GO ----------------
go_results <- list()
for (set in unique(gene_lists$set)) {
  genes <- list_genes(set)
  if (length(genes) < 2) {
    message("Lista '", set, "' com menos de 2 genes no universo: ignorada")
    next
  }
  for (ont in c("BP", "MF", "CC")) {
    res <- tryCatch(
      enrichGO(gene = genes, universe = universe, OrgDb = org.Hs.eg.db, keyType = key_type,
               ont = ont, pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
               minGSSize = min_size, maxGSSize = max_size, readable = FALSE),
      error = function(e) {
        message("enrichGO falhou (", set, ", ", ont, "): ", conditionMessage(e))
        NULL
      })
    out <- format_result(res, set, paste0("GO_", ont))
    if (!is.null(out)) go_results[[length(go_results) + 1]] <- out
  }
}
go_table <- if (length(go_results)) do.call(rbind, go_results) else empty_result
write.table(go_table, "GO.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
message("GO: ", nrow(go_table), " termos significativos")

# ---------------- KEGG ----------------
kegg_table <- empty_result
if (run_kegg) {
  kegg_failed <- FALSE
  ent <- if (key_type == "ENTREZID") {
    data.frame(ENTREZID = universe, stringsAsFactors = FALSE)
  } else {
    tryCatch(bitr(universe, fromType = key_type, toType = "ENTREZID", OrgDb = org.Hs.eg.db),
             error = function(e) {
               message("Conversão para ENTREZID falhou: ", conditionMessage(e))
               NULL
             })
  }
  if (!is.null(ent) && nrow(ent) > 0) {
    id_map <- setNames(ent[[key_type]], ent$ENTREZID)
    kegg_results <- list()
    for (set in unique(gene_lists$set)) {
      genes <- list_genes(set)
      ent_genes <- unique(ent$ENTREZID[ent[[key_type]] %in% genes])
      if (length(ent_genes) < 2) next
      res <- tryCatch(
        enrichKEGG(gene = ent_genes, universe = unique(ent$ENTREZID), organism = "hsa",
                   pvalueCutoff = 1, qvalueCutoff = 1,
                   minGSSize = min_size, maxGSSize = max_size),
        error = function(e) {
          message("enrichKEGG falhou (acesso à internet/KEGG?): ", conditionMessage(e))
          kegg_failed <<- TRUE
          NULL
        })
      if (kegg_failed) break
      out <- format_result(res, set, "KEGG", id_map)
      if (!is.null(out)) kegg_results[[length(kegg_results) + 1]] <- out
    }
    if (length(kegg_results)) kegg_table <- do.call(rbind, kegg_results)
  }
  message("KEGG: ", nrow(kegg_table), " vias significativas")
} else {
  message("KEGG desativado (--run_kegg false)")
}
write.table(kegg_table, "KEGG.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
