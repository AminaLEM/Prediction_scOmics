library(Seurat)
library(SeuratDisk)
set.seed(1234)


cohort_training = snakemake@input[[1]]
out1 = snakemake@output[[1]]

seu.filtered = readRDS(cohort_training)
# ch <-
#   NormalizeData(
#     object = ch,
#     normalization.method = "LogNormalize",
#     scale.factor = 10000
#   )
# RowsNA = rownames(ch)[rowSums(is.na(ch@assays$RNA@data)) > 0]
# '%!in%' <- function(x, y)
#   ! ('%in%'(x, y))
# RowsKEEP <- rownames(ch)[rownames(ch) %!in% RowsNA]
# ch <- subset(ch, features = RowsKEEP)
# ch$identif <-
#   paste(ch@meta.data$sampleID,
#         ch@meta.data$condition,
#         ch@meta.data$who_score,
#         ch@meta.data$batch,        
#         sep = "__")
# ch_avr = AverageExpression(
#   ch,
#   group.by = "identif",
#   verbose = FALSE,
#   return.seurat = TRUE,
#   slot = 'data'
# )
# ch_avr@meta.data[c('sampleID', 'condition', 'who_score','batch')] = t(data.frame(strsplit(colnames(ch_avr), "__")))
# AB.WT.index <- grep(pattern = "^AB-", x = rownames(ch_avr), value = FALSE) # Select row indices and not ERCC names 
# ch_avr =ch_avr[-AB.WT.index,]
# SaveH5Seurat(ch_avr, filename = out1, overwrite = TRUE)


DefaultAssay(seu.filtered)
seu.filtered$identif <-
  paste(seu.filtered@meta.data$sampleID,
        seu.filtered@meta.data$condition,
        seu.filtered@meta.data$who_score,
        seu.filtered@meta.data$batch,        
        sep = "__")

ch_avr <- AggregateExpression(seu.filtered, 
                           group.by = c("identif"),
                           assays = 'RNA',
                           slot = "counts",
                           return.seurat = TRUE)

ch_avr@meta.data[c('sampleID', 'condition', 'who_score','batch')] = t(data.frame(strsplit(colnames(ch_avr), "__")))
AB.WT.index <- grep(pattern = "^AB-", x = rownames(ch_avr), value = FALSE) # Select row indices and not ERCC names 
ch_avr =ch_avr[-AB.WT.index,]
SaveH5Seurat(ch_avr, filename = out1, overwrite = TRUE)
