library(Seurat)
library(SingleR)
library(celldex)
set.seed(0)
#  Load  QC data
cohort_training = snakemake@input[[1]]
out = snakemake@output[[1]]
pseud_mtx_all = NULL
ch = readRDS(cohort_training)
# Devide into batches
ch.batches <- SplitObject(ch, split.by =  "batch")
# Normalize per batch
ch.batches <-
  lapply(X = ch.batches, FUN = NormalizeData, normalization.method = "LogNormalize",scale.factor = 10000,verbose = FALSE)
'%!in%' <- function(x, y)   ! ('%in%'(x, y))
# Load monaco dataset for annotation
monaco.se <- MonacoImmuneData()
# Annotate per batch than merge cell composition matrices
for (ch in ch.batches) {
input <- GetAssayData(object = ch, slot = "data", assay = "RNA")
monaco <- SingleR(test = input, 
                               method="single",
                               fine.tune=FALSE,
                               ref = monaco.se, 
                               labels = monaco.se$label.fine)

ch$monaco.labels <- monaco$labels
ch$sampleID=factor(ch$sampleID) 
abundances <- table(ch$monaco.labels, ch$sampleID)
abundances <- unclass(abundances)
extra.info <-
  ch@meta.data[match(colnames(abundances), ch$sampleID), ]
rownames(extra.info) = colnames(abundances)
abundances = as.data.frame(t(abundances))
outer = unique(monaco.se$label.fine)[unique(monaco.se$label.fine) %!in% colnames(abundances)]
abundances[outer] = 0
abundances = abundances[, sort(names(abundances))]
abundances['condition'] = extra.info$condition
pseud_mtx_all = rbind(pseud_mtx_all, abundances)
    }
write.csv(pseud_mtx_all, paste0(out))