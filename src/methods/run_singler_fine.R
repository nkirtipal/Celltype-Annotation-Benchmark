suppressPackageStartupMessages({
  library(Matrix); library(SingleCellExperiment)
  library(scuttle); library(SingleR); library(celldex)
})

ROOT <- "/media/nikhil/My_Book/celltype_bench"
IN   <- file.path(ROOT, "data/processed/for_r")

counts   <- readMM(file.path(IN, "counts.mtx"))
rownames(counts) <- readLines(file.path(IN, "genes.tsv"))
colnames(counts) <- readLines(file.path(IN, "barcodes.tsv"))

sce <- SingleCellExperiment(assays = list(counts = counts))
sce <- logNormCounts(sce)

ref  <- celldex::HumanPrimaryCellAtlasData()
pred <- SingleR(test = sce, ref = ref, labels = ref$label.fine)

out <- data.frame(
  barcode    = rownames(pred),
  label      = pred$pruned.labels,
  confidence = apply(pred$scores, 1, max)
)
out$label[is.na(out$label)] <- "unassigned"

write.csv(out, file.path(ROOT, "results/predictions_singler_hpca_fine_pbmc3k.csv"),
          row.names = FALSE)

cat("n labels:", length(unique(out$label)), "\n")
print(sort(table(out$label), decreasing = TRUE)[1:25])
