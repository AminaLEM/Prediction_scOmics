library(Seurat)
library(SingleR)
library(celldex)
library(patchwork)
set.seed(0)
#  Load  QC data
cohort_training = snakemake@input[[1]]
out = snakemake@output[[1]]
pseud_mtx_all = NULL
ch = readRDS(cohort_training)
# Devide into batches
ch.batches <- SplitObject(ch, split.by =  "batch")
# Normalize per batch
ch.batches <-
  lapply(
    X = ch.batches,
    FUN = NormalizeData,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )
'%!in%' <- function(x, y)
  ! ('%in%'(x, y))
# Load monaco dataset for annotation
monaco.se <- MonacoImmuneData()
# Annotate per batch than merge cell composition matrices
for (ch in ch.batches) {
  input <- GetAssayData(object = ch,
                        slot = "data",
                        assay = "RNA")
  monaco <- SingleR(
    test = input,
    method = "single",
    fine.tune = FALSE,
    ref = monaco.se,
    labels = monaco.se$label.fine
  )
  
  ch$monaco.labels <- monaco$labels
  
  ch$sampleID = factor(ch$sampleID)
  ch$monaco.labels = factor(ch$monaco.labels)
  abundances <- table(ch$monaco.labels, ch$sampleID)
  abundances <- unclass(abundances)
  extra.info <-
    ch@meta.data[match(colnames(abundances), ch$sampleID),]
  rownames(extra.info) = colnames(abundances)
  abundances = as.data.frame(t(abundances))
  #Check for non expressed cells to add a zeros columns  
  outer = unique(monaco.se$label.fine)[unique(monaco.se$label.fine) %!in% colnames(abundances)]
  abundances[outer] = 0
  abundances = abundances[, sort(names(abundances))]
  abundances['condition1'] = extra.info$condition
  pseud_mtx_all = rbind(pseud_mtx_all, abundances)
}
# Merge cells into one celltype B cells, CD4 T cells, CD8 T cells, gd T cells, and change progenetors to Platelet cells
pseud_mtx_all['B cells'] = rowSums(pseud_mtx_all[, c("Naive B cells",
                                                     'Non-switched memory B cells',
                                                     'Switched memory B cells')], na.rm = TRUE)

pseud_mtx_all['CD4 T cells'] = rowSums(pseud_mtx_all[, c(
  "Th1/Th17 cells" ,
  "Th17 cells" ,
  "Th2 cells",
  "Follicular helper T cells",
  "Th1 cells",
  "Terminal effector CD4 T cells",
  "Naive CD4 T cells"
)], na.rm = TRUE)

pseud_mtx_all["CD8 T cells"] = rowSums(pseud_mtx_all[, c(
  "Terminal effector CD8 T cells",
  "Effector memory CD8 T cells",
  "Central memory CD8 T cells",
  "Naive CD8 T cells"
)], na.rm = TRUE)

pseud_mtx_all['gd T cells'] = rowSums(pseud_mtx_all[, c("Vd2 gd T cells", "Non-Vd2 gd T cells")], na.rm =
                                        TRUE)
pseud_mtx_all['Platelet cells'] = pseud_mtx_all[, "Progenitor cells"]


pseud_mtx_all['condition'] = pseud_mtx_all['condition1']
pseud_mtx_all = pseud_mtx_all[, !(
  names(pseud_mtx_all) %in% c(
    "Naive CD8 T cells",
    "Naive CD4 T cells",
    "Naive B cells",
    'Non-switched memory B cells',
    'Switched memory B cells',
    "Th1/Th17 cells" ,
    "Th17 cells" ,
    "Th2 cells",
    "Follicular helper T cells",
    "Th1 cells",
    "Terminal effector CD4 T cells",
    "Terminal effector CD8 T cells",
    "Effector memory CD8 T cells",
    "Central memory CD8 T cells",
    "Vd2 gd T cells",
    "Non-Vd2 gd T cells",
    'condition1',
    "Low-density basophils",
    "Progenitor cells",
    "Exhausted B cells"
  )
)]
write.csv(pseud_mtx_all, paste0(out))