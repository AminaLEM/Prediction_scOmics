library(SeuratDisk)
library(tidyverse)
set.seed(0)
# pseudo-bulk of the selected genes
inp = snakemake@input[[1]]
fc_file = snakemake@input[[2]]
out = snakemake@output[[1]]
ch_avr = LoadH5Seurat(paste0(inp,"_norm.h5Seurat"))
top = read.csv(fc_file, row.names = 1)[,1]
'%!in%' <- function(x, y)
  ! ('%in%'(x, y))
outer= top[top %!in% rownames(ch_avr)]
inner= top[top %in% rownames(ch_avr)]
pseud_mtx = as.data.frame(t(as.matrix(ch_avr@assays$RNA@data[inner, ])))
pseud_mtx[outer]=0
pseud_mtx = pseud_mtx[, sort(names(pseud_mtx))]
rownames(pseud_mtx) = ch_avr$sampleID
pseud_mtx['condition'] = ch_avr$condition
pseud_mtx['who_score'] = ch_avr$who_score
write.csv(pseud_mtx, out)