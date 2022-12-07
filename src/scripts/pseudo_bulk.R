# pseudo bulk fo data(for the downstream analysis) and count slots(to perform DE analysis)
library(Seurat)
library(SeuratDisk)
set.seed(1234)


cohort_training = snakemake@input[[1]]
out1 = snakemake@output[[1]]

seu.filtered = readRDS(cohort_training)
seu.filtered <-
  NormalizeData(
    object = seu.filtered,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )
DefaultAssay(seu.filtered)
seu.filtered$identif <-
  paste(
    seu.filtered@meta.data$sampleID,
    seu.filtered@meta.data$condition,
    seu.filtered@meta.data$who_score,
    seu.filtered@meta.data$batch,
    sep = "__"
  )

ch_avr <- AggregateExpression(
  seu.filtered,
  group.by = c("identif"),
  assays = 'RNA',
  slot = "counts",
  return.seurat = TRUE
)
ch_avr_nor <- AggregateExpression(
  seu.filtered,
  group.by = c("identif"),
  assays = 'RNA',
  slot = "data",
  return.seurat = TRUE
)

ch_avr@meta.data[c('sampleID', 'condition', 'who_score', 'batch')] = t(data.frame(strsplit(colnames(ch_avr), "__")))
ch_avr_nor@meta.data[c('sampleID', 'condition', 'who_score', 'batch')] = t(data.frame(strsplit(colnames(ch_avr_nor), "__")))

AB.WT.index <-
  grep(pattern = "^AB-",
       x = rownames(ch_avr),
       value = FALSE) 
ch_avr = ch_avr[-AB.WT.index, ]
SaveH5Seurat(ch_avr, filename = out1, overwrite = TRUE)
SaveH5Seurat(ch_avr_nor,
             filename = paste0(out1, "_norm.h5Seurat"),
             overwrite = TRUE)
