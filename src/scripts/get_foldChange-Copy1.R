library(Seurat)
library(SeuratDisk)
library(DESeq2)
library(tidyverse)
set.seed(0)
inp = snakemake@input[[1]]
out = snakemake@output[[1]]

ch = LoadH5Seurat(inp)
# Idents(ch) = 'condition'
# markers <-
#   FindMarkers(ch,
#               ident.1 = "Severe",
#               ident.2 = "Mild",
#               verbose = FALSE, 
#               latent.vars = "batch")


cts <- ch$RNA@counts



# perform DESeq2 --------
# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = ch@meta.data,
                              design = ~ batch+condition)

# filter
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]

# run DESeq2
dds <- DESeq(dds)

# Check the coefficients for the comparison
resultsNames(dds)

# Generate results object
markers <- results(dds, name = "condition_Severe_vs_Mild")


write.csv(as.data.frame(markers), out)