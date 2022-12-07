library(Seurat)
library(SeuratDisk)
library(DESeq2)
library(tidyverse)
library(edgeR)

set.seed(0)
# Load aggregated counts to sample level
inp <- snakemake@input[[1]]
out1 <- snakemake@output[[1]]
out2 <- snakemake@output[[2]]
out3 <- snakemake@output[[3]]
nb <- snakemake@params[[1]]

ch <- LoadH5Seurat(inp)
cts <- ch$RNA@counts
# Run edgeR --------
y <- DGEList(cts, samples = ch@meta.data)
keep <- filterByExpr(y, group = ch@meta.data$condition)
y <- y[keep, ]
y <- calcNormFactors(y)
design <- model.matrix( ~ 0 + batch + factor(condition), y$samples)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust = TRUE)
resedger <- glmQLFTest(fit, coef = ncol(design))
#Filter genes with |logFC|>2 and FDR<0.05
res3 = topTags(resedger, n = sum(abs(resedger$table$logFC) >= 2), sort.by =
                 "logFC")$table
res3 = res3[res3$FDR < 0.05, ]

# perform DESeq2 --------
# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = ch@meta.data,
  design = ~ batch + condition
)
# filter
dds <- dds[keep, ]
# run DESeq2
dds <- DESeq(dds)
# Generate results object
markers <- results(dds, name = "condition_Severe_vs_Mild")
res2 = as.data.frame(markers)
#Filter genes with |logFC|>2 and FDR<0.05
res2 = res2[abs(res2$log2FoldChange) >= 2 & res2$padj < 0.05, ]

# y <- DGEList(cts,samples=ch@meta.data)
# y <- calcNormFactors(y)

# design <- model.matrix(~0 + condition+batch, y$samples)

# cutoff <- 1
# drop <- which(apply(cpm(y), 1, max) < cutoff)
# d <- y[-drop,]
# dim(d) # number of genes left

# yy <- voom(d, design, plot = T)
# fit <- lmFit(yy, design)
# head(coef(fit))
# contr <- makeContrasts(conditionSevere-conditionMild, levels = colnames(coef(fit)))
# tmp <- contrasts.fit(fit, contr)
# tmp <- eBayes(tmp)
# top.table <- topTable(tmp, sort.by = "P", n = Inf)
# res=top.table[abs(top.table$logFC)>2.5,]
# res= res[res$adj.P.Val<0.05,]

# select intersection between TOPnb (from config file) of filtrerd genes from DESeq and edgeR
if (nb < dim(res2)[1]) {
  top_desq =  res2 %>% top_n(n = nb, wt = -abs(padj))
} else
{
  top_desq = res2
}
if (nb < dim(res3)[1]) {
  top_edgeR =  res3 %>% top_n(n = nb, wt = -abs(FDR))
} else
{
  top_edgeR = res3
}
selected = rownames(top_desq)[rownames(top_desq) %in% rownames(top_edgeR)]

#Save markers, counts and ...
write.csv(as.data.frame(markers), out1)
write.csv(as.data.frame(resedger$table), out2)
write.csv(as.data.frame(selected), out3)