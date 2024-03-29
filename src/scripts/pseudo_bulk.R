library(Seurat)
library(SeuratDisk)
cohort_training = snakemake@input[[1]]
out1 = snakemake@output[[1]]

ch = readRDS(cohort_training)
ch <-
  NormalizeData(
    object = ch,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )
RowsNA = rownames(ch)[rowSums(is.na(ch@assays$RNA@data)) > 0]
'%!in%' <- function(x, y)
  ! ('%in%'(x, y))
RowsKEEP <- rownames(ch)[rownames(ch) %!in% RowsNA]
ch <- subset(ch, features = RowsKEEP)
ch$identif <-
  paste(ch@meta.data$sampleID,
        ch@meta.data$condition,
        ch@meta.data$who_score,
        ch@meta.data$batch,        
        sep = "__")
ch_avr = AverageExpression(
  ch,
  group.by = "identif",
  verbose = FALSE,
  return.seurat = TRUE,
  slot = 'data'
)
ch_avr@meta.data[c('sampleID', 'condition', 'who_score','batch')] = t(data.frame(strsplit(colnames(ch_avr), "__")))
AB.WT.index <- grep(pattern = "^AB-", x = rownames(ch_avr), value = FALSE) # Select row indices and not ERCC names 
ch_avr =ch_avr[-AB.WT.index,]
SaveH5Seurat(ch_avr, filename = out1, overwrite = TRUE)
