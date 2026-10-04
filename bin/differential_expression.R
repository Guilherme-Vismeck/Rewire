#!/usr/bin/env Rscript
# Expressão diferencial com limma (DIFFERENTIAL_EXPRESSION)
suppressPackageStartupMessages(library(limma))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  i <- which(args == paste0("--", name))
  if (length(i) == 0) return(default)
  args[i[1] + 1]
}

expression_file <- get_arg("expression")
metadata_file   <- get_arg("metadata")
group_col       <- get_arg("group_col", "condition")
sample_col      <- get_arg("sample_column", "sample")
group_a         <- get_arg("group_a")
group_b         <- get_arg("group_b")
de_fdr          <- as.numeric(get_arg("de_fdr", "0.05"))
min_logfc       <- as.numeric(get_arg("min_logfc", "0"))
output          <- get_arg("output", "DE_results.tsv")

if (is.null(expression_file) || is.null(metadata_file) || is.null(group_a) || is.null(group_b)) {
  stop("Uso: differential_expression.R --expression X --metadata Y --group_a A --group_b B")
}

expr <- read.delim(expression_file, row.names = 1, check.names = FALSE)
meta <- read.delim(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)

keep <- meta[[group_col]] %in% c(group_a, group_b) & meta[[sample_col]] %in% colnames(expr)
meta <- meta[keep, , drop = FALSE]
expr <- as.matrix(expr[, meta[[sample_col]], drop = FALSE])

group <- factor(meta[[group_col]], levels = c(group_a, group_b))
n_a <- sum(group == group_a)
n_b <- sum(group == group_b)
if (n_a < 2 || n_b < 2) stop("Cada grupo precisa de pelo menos 2 amostras para o limma")

design <- model.matrix(~ group)
fit <- eBayes(lmFit(expr, design))
res <- topTable(fit, coef = 2, number = Inf, sort.by = "none")

out <- data.frame(
  gene    = rownames(res),
  logFC   = res$logFC,
  AveExpr = res$AveExpr,
  t       = res$t,
  pvalue  = res$P.Value,
  FDR     = res$adj.P.Val,
  stringsAsFactors = FALSE
)
out$DEG <- !is.na(out$FDR) & out$FDR < de_fdr & abs(out$logFC) >= min_logfc
out <- out[order(out$pvalue), ]

write.table(out, output, sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("Genes testados: %d (%s: %d amostras, %s: %d amostras)\n",
            nrow(out), group_a, n_a, group_b, n_b))
cat(sprintf("DEGs (FDR < %s e |logFC| >= %s): %d\n", de_fdr, min_logfc, sum(out$DEG)))
