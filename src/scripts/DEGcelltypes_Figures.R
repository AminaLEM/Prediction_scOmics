library(ggplot2)
library(reshape2)
library(Seurat)
library(SeuratDisk)
library(DESeq2)
library(tidyverse)
library(edgeR)

set.seed(0)

path_to_inp= snakemake@params[[1]]
path_to_outFig= snakemake@params[[2]]


# Load aggregated counts to sample level
cohort <-
  LoadH5Seurat(
    paste0(path_to_inp,"/merged_training/pseudo_bulk.h5Seurat")
  )
cohort$celltype = gsub('/', '', cohort$celltype)
cohort$celltype = gsub('Progenitor cells', 'Platelet cells', cohort$celltype)
# plot CPM for 15 random samples
all_celltypes = unique(cohort$celltype)[!unique(cohort$celltype) %in% c("Exhausted B cells", "Low-density basophils")]
all_results = NULL
for (ct in all_celltypes) {
  print(ct)
  ch = subset(cohort, subset = celltype == ct)
  cts <- ch$RNA@counts
  cts_cpm = cpm(cts)
  num <- seq_len(10)
  result <- data.frame(score = num,
                       do.call(rbind, lapply(num, function(x)
                         colSums(cts_cpm >= x))))
  result <-
    melt(result ,  id.vars = 'score', variable.name = 'series')
  result['series'] = substring(result$series, 1, regexpr("__", result$series) - 1)
  
  if (is.null(all_results)) {
    samples = sample(unique(result$series), 15)
    result = result[result$series %in% samples, ]
    result$celltype = ct
    all_results = result
  } else{
    result = result[result$series %in% samples, ]
    result$celltype = ct
    all_results = rbind(all_results, result)
  }
}

pdf(
  file = paste0(path_to_outFig,"/cpm.pdf"
  ),
  # The directory you want to save the file in
  width = 12,
  # The width of the plot in inches
  height = 10
) # The height of the plot in inches

ggplot(all_results, aes(score, value)) + geom_line(aes(colour = series)) +
  scale_colour_viridis_d() +
  xlab("CPM") + ylab("# of Genes > CPM Cutoff") + labs(color = 'Samples') +
  facet_wrap( ~ celltype, nrow = 3) + theme_bw() + theme(legend.justification = c(1, 0),
                                                         legend.text = element_text(size = 7)) + guides(color = guide_legend(ncol =
                                                                                                                               2))
dev.off()

#PLOT logFC
for (ct in all_celltypes) {
  print(ct)
  if (ct %in% c('Classical monocytes','gd T cells')){dimension = c(-2.5, 6)
  hg=5
  }else{
    dimension = c(-2, 3)  
    hg=4
    }
  ch = subset(cohort, subset = celltype == ct)
  cts <- ch$RNA@counts
  # Run edgeR --------
  y <-
    DGEList(cts,
            samples = ch@meta.data,
            group = ch@meta.data$condition)
  keep <- filterByExpr(y, group = ch@meta.data$condition)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  
  pdf(file = paste0(
    paste0(
      path_to_outFig,"/edgeR_",
      ct
    ),
    '.pdf'
  ),
  # The directory you want to save the file in
  width = 10,
  # The width of the plot in inches
  height = 10) # The height of the plot in inches
  
  par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
  pch <- c(0, 17)
  group = factor(ch@meta.data$condition)
  colors <- as.numeric(group)
  plotMDS(y,
          col = colors,
          pch = pch[group],
          main = ct)
  legend(
    "topleft",
    col = c(1, 2),
    legend = levels(group),
    pch = pch,
    ncol = 2
  )
  plotMD(y, column = 1, main = ct)
  abline(
    h = 0,
    col = "red",
    lty = 2,
    lwd = 2
  )
  
  design <- model.matrix( ~ 0 + batch + factor(condition), y$samples)
  y <- estimateDisp(y, design, robust = TRUE)
  plotBCV(y, main = ct)
  
  fit <- glmQLFit(y, design, robust = TRUE)
  plotQLDisp(fit, main = ct)
  dev.off()
  
  summary(fit$df.prior)
  resedger <- glmQLFTest(fit, coef = ncol(design))
  write.csv(as.data.frame(resedger$table), paste0(paste0(path_to_outFig,'/edgeR'),ct))
  
  is.de <- decideTestsDGE(resedger, p.value = 0.05)
  #Filter genes with |logFC|>2 and FDR<0.05
  nbDE = sum(abs(resedger$table$logFC) >= 2)
  if (nbDE > 0) {
    res3 = as.data.frame(topTags(resedger, n = sum(abs(
      resedger$table$logFC
    ) >= 2), sort.by = "logFC"))
  }
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
  write.csv(res2, paste0(paste0(path_to_outFig,'/deseq'),ct))
  res2 = res2[abs(res2$log2FoldChange) >= 2 & res2$padj < 0.05, ]
  
  pdf(file = paste0(
    paste0(paste0(
      path_to_outFig,'/'),
      ct
    ),
    '.pdf'
  ),
  # The directory you want to save the file in
  width = 10,
  # The width of the plot in inches
  height = hg) # The height of the plot in inches
  
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
    ylim = dimension,
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
  abline(
    h = 0,
    col = "#6b6b6c",
    lty = 1,
    lwd = 3
  )
  DESeq2::plotMA(markers, ylim = dimension, main = 'DESeq2')
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
  mtext(ct,
        outer = TRUE,
        cex = 1.5,
        font = 2)
  dev.off()
}
