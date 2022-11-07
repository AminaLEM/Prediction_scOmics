library(Seurat)
library(SeuratDisk)
set.seed(0)

cohort = snakemake@input[[1]]
out = snakemake@output[[1]]

ch = LoadH5Seurat(cohort)
ch[["percent.mt"]] <- PercentageFeatureSet(ch, pattern = "ˆMT-")
tryCatch(
  exp = {
    ch[['percent.HBB']] = PercentageFeatureSet(ch, features = c('HBB'), assay = 'RNA')
  },
  error = function(e) {
    print('no HBB gene ...')
  },
  finally = {
    ch[['percent.HBB']] = 0
  }
)
tryCatch(
  exp = {
    ch[['percent.HBA']] = PercentageFeatureSet(ch, features = c('HBA'), assay = 'RNA')
  },
  error = function(e) {
    print('no HBA gene ...')
  },
  finally = {
    ch[['percent.HBA']] = 0
  }
)
ch <-
  subset(
    ch,
    subset = nFeature_RNA > 250 &
      nFeature_RNA < 5000 &
      percent.mt < 25 &
      nCount_RNA > 500 & percent.HBB < 25 & percent.HBA < 25
  )
AB.WT.index <- grep(pattern = "^AB-", x = rownames(ch), value = FALSE) 
ch =ch[-AB.WT.index,]

saveRDS(ch, out)
