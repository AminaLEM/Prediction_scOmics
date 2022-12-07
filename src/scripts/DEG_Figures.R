library(Seurat)
library(SeuratDisk)
library(DESeq2)
library(tidyverse)
library(edgeR)
library(reshape2)

library(pheatmap)
library(RColorBrewer)
set.seed(0)
path_to_inp= snakemake@params[[1]]
path_to_inpFig= snakemake@params[[2]]

inp <-
  paste0(path_to_inp,"/merged_training/pseudo_bulk.h5Seurat")
ch <- LoadH5Seurat(inp)


set.seed(0)
cts <- ch$RNA@counts
cts_cpm = cpm(cts)
num <- seq_len(10)
result <- data.frame(score = num,
                     do.call(rbind, lapply(num, function(x)
                       colSums(cts_cpm >= x))))
result <-
  melt(result ,  id.vars = 'score', variable.name = 'series')
result['series'] = substring(result$series, 1, regexpr("__", result$series) - 1)
samples = sample(unique(result$series), 15)
result = result[result$series %in% samples, ]

pdf(
  file = paste0(path_to_inpFig,"/cpm_all.pdf"),
  # The directory you want to save the file in
  width = 12,
  # The width of the plot in inches
  height = 10
) # The height of the plot in inches

ggplot(result, aes(score, value)) + geom_line(aes(colour = series)) +
  scale_colour_viridis_d() +
  xlab("CPM") + ylab("# of Genes > CPM Cutoff") + labs(color = 'Samples') +
  theme_bw() + theme(legend.justification = c(1, 0),
                     legend.text = element_text(size = 7)) + guides(color = guide_legend(ncol =
                                                                                           2))
dev.off()


#PLOT logFC
cts <- ch$RNA@counts
# Run edgeR --------
y <- DGEList(cts,
             samples = ch@meta.data,
             group = ch@meta.data$condition)
keep <- filterByExpr(y, group = ch@meta.data$condition)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

pdf(
  file = paste0(path_to_inpFig,"/edgeR_all.pdf"
  ),
  # The directory you want to save the file in
  width = 10,
  # The width of the plot in inches
  height = 10
) # The height of the plot in inches

par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
pch <- c(0, 17)
group = factor(ch@meta.data$condition)
colors <- as.numeric(group)
plotMDS(y, col = colors, pch = pch[group])
legend(
  "topleft",
  col = c(1, 2),
  legend = levels(group),
  pch = pch,
  ncol = 2
)
plotMD(y, column = 1)
abline(h = 0,
       col = "red",
       lty = 2,
       lwd = 2)

design <- model.matrix( ~ 0 + batch + factor(condition), y$samples)
y <- estimateDisp(y, design, robust = TRUE)
plotBCV(y)

fit <- glmQLFit(y, design, robust = TRUE)
plotQLDisp(fit)
dev.off()

summary(fit$df.prior)
resedger <- glmQLFTest(fit, coef = ncol(design))
is.de <- decideTestsDGE(resedger, p.value = 0.05)
#Filter genes with |logFC|>2 and FDR<0.05
nbDE = sum(abs(resedger$table$logFC) >= 2)
if (nbDE > 0) {
  res3 = as.data.frame(topTags(resedger, n = sum(abs(
    resedger$table$logFC
  ) >= 2), sort.by = "logFC"))
}
#$table

#tr <- glmTreat(fit, coef=ncol(design), lfc=log2(1.5))
#is.de <- decideTestsDGE(tr)
#plotMD(resedger, status=is.de, values=c(1,-1), col=c("red","blue"), legend="topright", main=ct)

# perform DESeq2 --------
# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = ch@meta.data,
  design = ~ batch + condition
)
# filter
#keep <- rowSums(counts(dds) >= 10) >= 10
dds <- dds[keep, ]
# run DESeq2
dds <- DESeq(dds)
# Generate results object
markers <-
  results(dds, name = "condition_Severe_vs_Mild", alpha = 0.05)
res2 = as.data.frame(markers)
res2 = res2[abs(res2$log2FoldChange) >= 2 & res2$padj < 0.05, ]

pdf(
  file = paste0(path_to_inpFig,"/DEG_all.pdf"),
  # The directory you want to save the file in
  width = 15,
  # The width of the plot in inches
  height = 10
) # The height of the plot in inches

par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))
plotMD(
  resedger,
  hl.cex = 0.3,
  bg.col = '#86868e',
  xlab = 'average log CPM',
  ylab = 'log fold change',
  status = is.de,
  values = c(1, -1),
  col = c("blue", "blue"),
  legend = FALSE,
  main = 'edgeR',
  ylim = c(-4, 6.5),
  pch = 20
)
if (nbDE > 0) {
  text(
    res3$logCPM,
    res3$logFC,
    labels = rownames(res3),
    col = "black",
    cex = 0.7,
    pos = 3,
    font = 2
  )
}
abline(h = 0,
       col = "#6b6b6c",
       lty = 1,
       lwd = 3)
DESeq2::plotMA(markers, ylim = c(-4, 6.5), main = 'DESeq2')
if (dim(res2)[1] > 0) {
  text(
    res2$baseMean,
    res2$log2FoldChange,
    labels = rownames(res2),
    col = "black",
    cex = 0.7,
    pos = 3,
    font = 2
  )
}
#mtext(ct, outer = TRUE, cex = 1.5,font=2)
dev.off()
# top_edgeR =  res3 %>% top_n(n = 10, wt = -abs(FDR))
# top_desq =  res2 %>% top_n(n = 10, wt = -abs(padj))

# selected = rownames(top_edgeR)[rownames(top_edgeR) %in% rownames(top_desq)]



# mydf = cpm(y, normalized.lib.sizes = FALSE)[selected, ]
# my_sample_col <- data.frame(condition = group, batch = ch$batch)
# row.names(my_sample_col) <- ch$sampleID
# my_sample_col = arrange(my_sample_col, condition)
# my_sample_col = my_sample_col[my_sample_col$batch == '1_trainset_2', ]
# colnames(mydf) = ch$sampleID

# mydf = mydf[, row.names(my_sample_col)]
# pheatmap(
#   mydf,
#   cluster_rows = F,
#   cluster_cols = F,
#   annotation_col = my_sample_col
# )


# mydf = cpm(y, normalized.lib.sizes = FALSE)[selected, ]
# my_sample_col <- data.frame(condition = group, batch = ch$batch)
# row.names(my_sample_col) <- ch$sampleID
# my_sample_col = arrange(my_sample_col, condition)
# my_sample_col = my_sample_col[my_sample_col$batch == '1_trainset_1', ]
# colnames(mydf) = ch$sampleID

# mydf = mydf[, row.names(my_sample_col)]
# pheatmap(
#   mydf[c('RETN', 'S100P', 'ANXA3'), ],
#   cluster_rows = F,
#   cluster_cols = F,
#   annotation_col = my_sample_col
# )
