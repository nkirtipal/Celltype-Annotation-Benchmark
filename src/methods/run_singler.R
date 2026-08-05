suppressPackageStartupMessages({
  library(Matrix); library(SingleCellExperiment)
  library(scuttle); library(SingleR); library(celldex)
})

ROOT <- "/media/nikhil/My_Book/celltype_bench"
IN   <- file.path(ROOT, "data/processed/for_r")

counts   <- readMM(file.path(IN, "counts.mtx"))
genes    <- readLines(file.path(IN, "genes.tsv"))
barcodes <- readLines(file.path(IN, "barcodes.tsv"))

rownames(counts) <- genes
colnames(counts) <- barcodes
cat("matrix:", dim(counts), "\n")

sce <- SingleCellExperiment(assays = list(counts = counts))
sce <- logNormCounts(sce)

ref  <- celldex::HumanPrimaryCellAtlasData()
cat("reference:", dim(ref), "\n")

pred <- SingleR(test = sce, ref = ref, labels = ref$label.main)

out <- data.frame(
  barcode    = rownames(pred),
  label      = pred$pruned.labels,
  confidence = apply(pred$scores, 1, max)
)
out$label[is.na(out$label)] <- "unassigned"

write.csv(out, file.path(ROOT, "results/predictions_singler_hpca_pbmc3k.csv"),
          row.names = FALSE)

cat("done\n")
print(table(out$label))
