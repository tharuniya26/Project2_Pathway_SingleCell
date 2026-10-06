# Differential expression analysis for beta-cell pseudobulk data
# Comparison: T2D vs ND
# Technical covariate: Chemistry

.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

library(DESeq2)

# Load exported data
counts <- read.csv(
  "data/processed/beta_counts.csv",
  row.names = 1,
  check.names = FALSE
)

design <- read.csv(
  "data/processed/beta_design.csv",
  row.names = 1,
  check.names = FALSE
)

# Keep only ND and T2D donors
design <- design[design$Condition %in% c("ND", "T2D"), ]

# Keep matching count columns in exactly the same order
counts <- counts[, rownames(design)]

# Confirm donor alignment
stopifnot(
  identical(colnames(counts), rownames(design))
)

# Set factor reference levels
design$Condition <- factor(
  design$Condition,
  levels = c("ND", "T2D")
)

design$Chemistry <- factor(design$Chemistry)

cat("Donors included:\n")
print(table(design$Condition))

cat("\nChemistry distribution:\n")
print(table(design$Chemistry, design$Condition))

# Filter low-expression genes
keep <- rowSums(counts >= 10) >= 5
counts_filtered <- counts[keep, ]

cat("\nGenes before filtering:", nrow(counts), "\n")
cat("Genes after filtering:", nrow(counts_filtered), "\n")

# DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts_filtered)),
  colData = design,
  design = ~ Chemistry + Condition
)

# Run DESeq2
dds <- DESeq(dds)

# T2D versus ND
res <- results(
  dds,
  contrast = c("Condition", "T2D", "ND"),
  alpha = 0.05
)

# Order by adjusted p-value
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[
  order(res_df$padj, na.last = TRUE),
]

# Summary
cat("\nDESeq2 result summary:\n")
print(summary(res))

cat(
  "\nSignificant genes (padj < 0.05):",
  sum(res_df$padj < 0.05, na.rm = TRUE),
  "\n"
)

cat(
  "Upregulated in T2D (padj < 0.05 and log2FC > 0):",
  sum(
    res_df$padj < 0.05 &
    res_df$log2FoldChange > 0,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Downregulated in T2D (padj < 0.05 and log2FC < 0):",
  sum(
    res_df$padj < 0.05 &
    res_df$log2FoldChange < 0,
    na.rm = TRUE
  ),
  "\n"
)

# Save full results
dir.create(
  "results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  res_df,
  "results/tables/beta_T2D_vs_ND_DESeq2.csv",
  row.names = FALSE
)

cat(
  "\nSaved results to:\n",
  "results/tables/beta_T2D_vs_ND_DESeq2.csv\n"
)