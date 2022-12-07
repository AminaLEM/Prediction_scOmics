library(Seurat)
library(future)
#\plan("multiprocess", workers = 20)
set.seed(0)

args = commandArgs(trailingOnly=TRUE)

ch = readRDS(args[1])

ch <-SCTransform(ch, vars.to.regress = c("percent.mt","nCount_RNA", "percent.HBB","percent.HBA"), verbose = FALSE,  return.only.var.genes = TRUE)
ch <- RunPCA(ch, verbose = TRUE)
ch <- RunUMAP(ch, min.dist = 0.15, dims = 1:50, verbose = F, reduction.name = "umap_mindist0.15", seed.use = 42)


Idents(ch)=ch$celltype
ch <- PrepSCTFindMarkers(ch)
markers <- FindAllMarkers(object = ch,assay = "SCT",
                                                 only.pos = TRUE,
                                                 min.pct = 0.2,
                                                 logfc.threshold = 0.2,
                                                 min.diff.pct = 0.1,
                                                 test.use = "wilcox"
)

write.csv(markers, paste0(args[1],'.csv'))
saveRDS(ch, args[1])