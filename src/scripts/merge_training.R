library(Seurat)
set.seed(0)
cohort = snakemake@input
n = length(cohort)
#Merge training sets and set as different batches
for (i in 1:n) {
  print(cohort[[i]])
  chrt = readRDS(cohort[[i]])
  chr_n = strsplit(cohort[[i]], "/")[[1]]
  chr_n = chr_n[length(chr_n)]
  chrt@meta.data['batch'] = paste0(chrt@meta.data$batch, "_trainset_", i)
  if (i==1) {
    cohort_merged = chrt
  }
  else {
    cohort_merged = merge(cohort_merged, chrt)
  }
}
saveRDS(cohort_merged, snakemake@output[[1]])