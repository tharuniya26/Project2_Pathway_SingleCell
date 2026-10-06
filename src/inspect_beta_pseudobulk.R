# Load the beta-cell pseudobulk object
x <- readRDS("data/processed/Rmarkdown/BetaPseudobulk.rds")

# Separate expression counts and donor metadata
counts <- x[["Counts"]]
design <- x[["Design"]]

# Show count matrix dimensions
cat("Counts dimensions:\n")
print(dim(counts))

# Show disease group counts
cat("\nDisease group counts:\n")
print(table(design$Condition))

# Show first donor IDs from the count matrix
cat("\nFirst count-matrix donor IDs:\n")
print(head(colnames(counts)))

# Show first donor IDs from the metadata
cat("\nFirst metadata donor IDs:\n")
print(head(rownames(design)))

# Check whether the count columns and metadata rows
# refer to the same donors in the same order
same_order <- identical(
  colnames(counts),
  rownames(design)
)

cat("\nSame donor order:", same_order, "\n")
# Export counts and donor metadata for Python analysis
write.csv(
  counts,
  "data/processed/beta_counts.csv",
  row.names = TRUE
)

write.csv(
  design,
  "data/processed/beta_design.csv",
  row.names = TRUE
)

cat("\nExport complete:\n")
cat("data/processed/beta_counts.csv\n")
cat("data/processed/beta_design.csv\n")