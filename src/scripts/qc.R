library(Seurat)
library(SeuratDisk)
set.seed(0)

cohort = snakemake@input[[1]]
out = snakemake@output[[1]]

ch = LoadH5Seurat(cohort)
ch[["percent.mt"]] <- PercentageFeatureSet(ch, pattern = "ˆMT-")
ch[["percent.rRNA28s"]] <-
  PercentageFeatureSet(ch, pattern = "RNA28S5")
ch[["percent.rRNA18s"]] <-
  PercentageFeatureSet(ch, pattern = "RNA18S5")
ch[["percent.HB"]] <- PercentageFeatureSet(ch, pattern = "^HB-")


ch <-
  subset(
    ch,
    subset = nFeature_RNA > 250 &
      nFeature_RNA / nCount_RNA < 0.75 &
      percent.mt < 20 &
      percent.rRNA28s < 20 &
      percent.rRNA18s < 20 &
      
      nCount_RNA > 700 & percent.HB < 20
      #nCount_RNA < 15000 & percent.HB < 20
  )
AB.WT.index <-
  grep(pattern = "^AB-",
       x = rownames(ch),
       value = FALSE)
ch = ch[-AB.WT.index, ]
saveRDS(ch, out)
