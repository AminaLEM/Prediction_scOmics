## ----setup, echo=FALSE, include=FALSE, message=FALSE--------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE, cache.lazy = FALSE)
knitr::opts_chunk$set(fig.width=15, fig.height=7)
knitr::opts_chunk$set(dev="CairoPNG")
knitr::opts_knit$set(root.dir = ".")
set.seed(123)


## ----libraries----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(Seurat)
library(future)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(SingleR)
library(celldex)
library(RColorBrewer)
library(SeuratDisk)
library(ggpubr)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bonn= readRDS('/prj/NUM_CODEX_PLUS/Amina/CellSubmission/new_runs/output_TOP5/bonn.h5Seurat/QC.rds')
berlin= readRDS('/prj/NUM_CODEX_PLUS/Amina/CellSubmission/new_runs/output_TOP5/berlin.h5Seurat/QC.rds')
kor= readRDS('/prj/NUM_CODEX_PLUS/Amina/CellSubmission/new_runs/output_TOP5/korean.h5Seurat/QC.rds')
stan= readRDS('/prj/NUM_CODEX_PLUS/Amina/CellSubmission/new_runs/output_TOP5/stanford.h5Seurat/QC.rds')


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

bonn <- NormalizeData(object = bonn, normalization.method = "LogNormalize", scale.factor = 10000)
berlin <- NormalizeData(object = berlin, normalization.method = "LogNormalize", scale.factor = 10000)
kor <- NormalizeData(object = kor, normalization.method = "LogNormalize", scale.factor = 10000)
stan <- NormalizeData(object = stan, normalization.method = "LogNormalize", scale.factor = 10000)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

bonn <- FindVariableFeatures(object = bonn, assay = "RNA", selection.method = "vst")
berlin <- FindVariableFeatures(object = berlin, assay = "RNA", selection.method = "vst")
kor <- FindVariableFeatures(object = kor, assay = "RNA", selection.method = "vst")
stan <- FindVariableFeatures(object = stan, assay = "RNA", selection.method = "vst", nfeatures = 3000)



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

top25_bonn <- head(VariableFeatures(bonn), 25)
top25_berlin <- head(VariableFeatures(berlin), 25)
top25_kor <- head(VariableFeatures(kor), 25)
top25_stan <- head(VariableFeatures(stan), 25)



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

plot1 <- VariableFeaturePlot(bonn)+ggtitle('Bonn cohort')
LabelPoints(plot = plot1, points = top25_bonn, repel = TRUE)

plot1 <- VariableFeaturePlot(berlin)+ggtitle('Berlin cohort')
LabelPoints(plot = plot1, points = top25_berlin, repel = TRUE)

plot1 <- VariableFeaturePlot(kor)+ggtitle('Korean cohort')
LabelPoints(plot = plot1, points = top25_kor, repel = TRUE)

plot1 <- VariableFeaturePlot(stan)+ggtitle('Stanford cohort')
LabelPoints(plot = plot1, points = top25_stan, repel = TRUE)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

bonn <- ScaleData(object = bonn, vars.to.regress = c("nCount_RNA"))
berlin <- ScaleData(object = berlin, vars.to.regress = c("nCount_RNA"))
kor <- ScaleData(object = kor, vars.to.regress = c("nCount_RNA"))
stan <- ScaleData(object = stan, vars.to.regress = c("nCount_RNA"))



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

bonn <- RunPCA(object = bonn, features = VariableFeatures(object = bonn),  verbose = FALSE)
berlin <- RunPCA(object = berlin, features = VariableFeatures(object = berlin),  verbose = FALSE)
kor <- RunPCA(object = kor, features = VariableFeatures(object = kor),  verbose = FALSE)
stan <- RunPCA(object = stan, features = VariableFeatures(object = stan),  verbose = FALSE)



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
DimHeatmap(bonn, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(berlin, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(kor, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(stan, dims = 1, cells = 500, balanced = TRUE)



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

ElbowPlot(bonn,ndims = 30)+ggtitle('bonn cohort')
ElbowPlot(berlin,ndims = 30)+ggtitle('berlin cohort')
ElbowPlot(kor,ndims = 30)+ggtitle('korean cohort')
ElbowPlot(stan,ndims = 50)+ggtitle('Stanford cohort')



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

bonn <- RunUMAP(bonn, reduction.use = "pca",  dims = 1:20, seed.use = 42)
berlin <- RunUMAP(berlin, reduction.use = "pca",  dims = 1:20, seed.use = 42)
kor <- RunUMAP(kor, reduction.use = "pca",  dims = 1:20, seed.use = 42)
stan <- RunUMAP(stan, reduction.use = "pca",  dims = 1:50, seed.use = 42)



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
DimPlot(object = bonn, 
        reduction = 'umap', 
        label = FALSE, 
        group.by = "sampleID"
        )
DimPlot(object = berlin, 
        reduction = 'umap', 
        label = FALSE, 
        group.by = "sampleID"
        )
DimPlot(object = kor, 
        reduction = 'umap', 
        label = FALSE, 
        group.by = "sampleID"
        )
DimPlot(object = stan, 
        reduction = 'umap', 
        label = FALSE, 
        group.by = "sampleID"
        )


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
monaco.se <- MonacoImmuneData()
#monaco.se <- NovershternHematopoieticData()


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# perform singleR classification
bonn$monaco.labels.f <- SingleR(test = GetAssayData(object = bonn, slot = "data", assay = "RNA"), 
                               method="single",
                               fine.tune=FALSE,
                               ref = monaco.se, 
                               labels = monaco.se$label.fine)$labels

berlin$monaco.labels.f <- SingleR(test = GetAssayData(object = berlin, slot = "data", assay = "RNA"), 
                               method="single",
                               fine.tune=FALSE,
                               ref = monaco.se, 
                               labels = monaco.se$label.fine)$labels
kor$monaco.labels.f <- SingleR(test = GetAssayData(object = kor, slot = "data", assay = "RNA"), 
                               method="single",
                               fine.tune=FALSE,
                               ref = monaco.se, 
                               labels = monaco.se$label.fine)$labels
stan$monaco.labels.f <- SingleR(test = GetAssayData(object = stan, slot = "data", assay = "RNA"), 
                               method="single",
                               fine.tune=FALSE,
                               ref = monaco.se, 
                               labels = monaco.se$label.fine)$labels


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bonn_ql <- SingleR(test = GetAssayData(object = bonn, slot = "data", assay = "RNA"), 
                               method="single",
                               fine.tune=FALSE,
                               ref = monaco.se, 
                               labels = monaco.se$label.fine)
p1=plotScoreHeatmap(bonn_ql)

berlin_ql <- SingleR(test = GetAssayData(object = berlin, slot = "data", assay = "RNA"), 
                               method="single",
                               fine.tune=FALSE,
                               ref = monaco.se, 
                               labels = monaco.se$label.fine)
p2=plotScoreHeatmap(berlin_ql)

kor_ql <- SingleR(test = GetAssayData(object = kor, slot = "data", assay = "RNA"), 
                               method="single",
                               fine.tune=FALSE,
                               ref = monaco.se, 
                               labels = monaco.se$label.fine)
p3=plotScoreHeatmap(kor_ql)

stan_ql<- SingleR(test = GetAssayData(object = stan, slot = "data", assay = "RNA"), 
                               method="single",
                               fine.tune=FALSE,
                               ref = monaco.se, 
                               labels = monaco.se$label.fine)
p4=plotScoreHeatmap(stan_ql)

png(file = paste0("~/new_NUM_CODEX_PLUS/Amina/CellSubmission/new_runs/figures/singlerRscores_bonn.pdf"),   # The directory you want to save the file in
    width = 15, # The width of the plot in inches
    height = 5) # The height of the plot in inches
p1
dev.off()

pdf(file = paste0("~/new_NUM_CODEX_PLUS/Amina/CellSubmission/new_runs/figures/singlerRscores_berlin.pdf"),   # The directory you want to save the file in
    width = 15, # The width of the plot in inches
    height = 5) # The height of the plot in inches
p2
dev.off()
pdf(file = paste0("~/new_NUM_CODEX_PLUS/Amina/CellSubmission/new_runs/figures/singlerRscores_kor.pdf"),   # The directory you want to save the file in
    width = 15, # The width of the plot in inches
    height = 5) # The height of the plot in inches
p3
dev.off()
pdf(file = paste0("~/new_NUM_CODEX_PLUS/Amina/CellSubmission/new_runs/figures/singlerRscores_stan.pdf"),   # The directory you want to save the file in
    width = 15, # The width of the plot in inches
    height = 5) # The height of the plot in inches
p4
dev.off()




## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
annotate <- function(bonn){
bonn$celltype = bonn$monaco.labels.f



bonn$celltype[bonn$celltype %in% c("Naive B cells",'Non-switched memory B cells','Switched memory B cells')] = 'B cells'
bonn$celltype[bonn$celltype %in% c("Th1/Th17 cells" ,"Th17 cells" ,"Th2 cells","Follicular helper T cells","Th1 cells","Terminal effector CD4 T cells","Naive CD4 T cells")] = 'CD4 T cells'
bonn$celltype[bonn$celltype %in% c("Terminal effector CD8 T cells","Effector memory CD8 T cells","Central memory CD8 T cells","Naive CD8 T cells")] = "CD8 T cells"

bonn$celltype[bonn$celltype %in% c("Vd2 gd T cells","Non-Vd2 gd T cells")] = 'gd T cells'
bonn$celltype[bonn$celltype == "Progenitor cells"] = 'Platelet cells'
bonn$celltype[bonn$celltype %in% c("Low-density basophils","Exhausted B cells")] = 'mixed'


identities = c('Classical monocytes',"Non classical monocytes","Intermediate monocytes", "Plasmacytoid dendritic cells","Myeloid dendritic cells", "B cells", "Plasmablasts", "CD8 T cells" , "CD4 T cells", "T regulatory cells", "MAIT cells", "gd T cells","Natural killer cells","Low-density neutrophils" ,"Platelet cells",'mixed' )
#identities = c('Classical monocytes',"Non classical monocytes","Intermediate monocytes", "Plasmacytoid dendritic cells","Myeloid dendritic cells", "Naive B cells","Exhausted B cells" , "memory B cells" ,"Plasmablasts", "CD8 T cells" , "CD4 T cells", "MAIT cells", "gd T cells","Natural killer cells","Low-density neutrophils","Low-density basophils" ,"Progenitor cells" )

#Idents(bonn)= bonn$celltype

Idents(bonn)=factor(bonn$celltype, levels = identities)
#DotPlot(bonn, features = c('LRRN3', 'CCR7', 'NPM1', 'SELL', 'PASK', 'IL7R', 'KLRB1', 'PRF1','TNFSF13B','GZMK','CCL5','CCL4','GZMH','GZMA','GNLY','NKG7','CST7','ITGB1','FOXP3','IL2RA','CD4','CD8A','CD8B','TRAC','CD3G','CD3D','CD3E','CCR5','TBX21','TGFB1'), cols = 'RdBu')+RotatedAxis()+ggtitle('stanford')
return(bonn)
}
bonn =annotate (bonn)
stan =annotate (stan)
kor =annotate (kor)
berlin =annotate (berlin)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
findmarker <- function(ch,dataset){
Idents(ch)=ch$celltype
markers <- FindAllMarkers(object = ch
                                                 only.pos = TRUE,
                                                 min.pct = 0.2,
                                                 logfc.threshold = 0.2,
                                                 min.diff.pct = 0.1,
                                                 test.use = "wilcox"
)

write.csv(markers, paste0(paste0('~/new_NUM_CODEX_PLUS/Amina/CellSubmission/new_runs/figures/markers_',dataset),'.csv'))
return(markers)
}
bonn_markers = findmarker(bonn,'bonn')
berlin_markers = findmarker(berlin,'berlin')
kor_markers = findmarker(kor,'korean')
stan_markers = findmarker(stan,'stanford')

bonn_markers = bonn_markers %>% filter(abs(avg_log2FC) > 1, p_val_adj <0.05)
berlin_markers = berlin_markers %>% filter(abs(avg_log2FC) > 1, p_val_adj <0.05)
kor_markers = kor_markers %>% filter(abs(avg_log2FC) > 1, p_val_adj <0.05)
stan_markers = stan_markers %>% filter(abs(avg_log2FC) > 1, p_val_adj <0.05)

bonn_markers <- bonn_markers %>% group_by(cluster) %>% top_n(n = 25, wt = abs(avg_log2FC))
berlin_markers <- berlin_markers %>% group_by(cluster) %>% top_n(n = 25, wt = abs(avg_log2FC))
kor_markers <- kor_markers %>% group_by(cluster) %>% top_n(n = 25, wt = abs(avg_log2FC))
stan_markers <- stan_markers %>% group_by(cluster) %>% top_n(n = 25, wt = abs(avg_log2FC))

merged = merge(
  x=bonn_markers,
  y=berlin_markers,
  by=c("cluster","gene"),

)
merged = merge(
  x=merged,
  y=kor_markers,
  by=c("cluster","gene"),

)
merged = merge(
  x=merged,
  y=stan_markers,
  by=c("cluster","gene"),
  suffixes = c(".merged",".stan")
)
#merged=merged %>% arrange(factor(cluster, levels = identities))
top <- merged %>% group_by(cluster) %>% top_n(n = 5, wt = abs(avg_log2FC.x))
top


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
identities = c('Classical monocytes',"Non classical monocytes","Intermediate monocytes", "Plasmacytoid dendritic cells","Myeloid dendritic cells", "B cells", "Plasmablasts", "CD8 T cells" , "CD4 T cells", "T regulatory cells", "MAIT cells", "gd T cells","Natural killer cells","Low-density neutrophils" ,"Platelet cells",'mixed' )
# 
# annotate <- function(bonn){
# Idents(bonn)=factor(bonn$celltype, levels = identities)
# return(bonn)
# }
# bonn =annotate (bonn)
# stan =annotate (stan)
# kor =annotate (kor)
# berlin =annotate (berlin)


bonn_markers = read.csv( '/prj/NUM_CODEX_PLUS/Amina/CellSubmission/code_v1/data/bonn_.rds.csv')
berlin_markers = read.csv( '/prj/NUM_CODEX_PLUS/Amina/CellSubmission/code_v1/data/berlin_.rds.csv')
kor_markers = read.csv( '/prj/NUM_CODEX_PLUS/Amina/CellSubmission/code_v1/data/kor_.rds.csv')
stan_markers = read.csv( '/prj/NUM_CODEX_PLUS/Amina/CellSubmission/code_v1/data/stan_v0.rds.csv')

bonn_markers = bonn_markers %>% filter(abs(avg_log2FC) > 1, p_val_adj <0.05)
berlin_markers = berlin_markers %>% filter(abs(avg_log2FC) > 1, p_val_adj <0.05)
kor_markers = kor_markers %>% filter(abs(avg_log2FC) > 1, p_val_adj <0.05)
stan_markers = stan_markers %>% filter(abs(avg_log2FC) > 1, p_val_adj <0.05)

bonn_markers <- bonn_markers %>% group_by(cluster) %>% top_n(n = 5, wt = abs(avg_log2FC))%>%  arrange(factor(cluster, levels = identities))
berlin_markers <- berlin_markers %>% group_by(cluster) %>% top_n(n = 5, wt = abs(avg_log2FC))%>%  arrange(factor(cluster, levels = identities))
kor_markers <- kor_markers %>% group_by(cluster) %>% top_n(n = 5, wt = abs(avg_log2FC))%>%  arrange(factor(cluster, levels = identities))
stan_markers <- stan_markers %>% group_by(cluster) %>% top_n(n = 5, wt = abs(avg_log2FC))%>%  arrange(factor(cluster, levels = identities))



p1 = DotPlot(bonn , features = unique(bonn_markers$gene)[unique(bonn_markers$gene)!='FCER1A'], cols = 'RdBu')+RotatedAxis()+theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"),        axis.text.x=element_text(angle=90, hjust=1, size=10),axis.title=element_blank()) +NoLegend()
p2= DotPlot(berlin , features = unique(berlin_markers$gene), cols = 'RdBu')+RotatedAxis()+theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"),        axis.text.x=element_text(angle=90, hjust=1, size=10),axis.title=element_blank())+NoLegend()
p4=DotPlot(kor , features = unique(kor_markers$gene), cols = 'RdBu')+RotatedAxis()+theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"),        axis.text.x=element_text(angle=90, hjust=1, size=10),axis.title=element_blank())
p3=DotPlot(stan , features = unique(stan_markers$gene), cols = 'RdBu')+RotatedAxis()+theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"),        axis.text.x=element_text(angle=90, hjust=1, size=10),axis.title=element_blank())+NoLegend()

figure <- ggarrange(p1,p2,p3,p4,
                    font.label = list(size = 16, color = "black", face = "bold", family = NULL),
                    labels = c("Bonn", "Berlin","Stanford",'Korea'),vjust= 1.2,
                    ncol = 1, nrow = 4)

pdf(file = paste0("~/new_NUM_CODEX_PLUS/Amina/CellSubmission/new_runs/figures/Topmarkergenes.pdf"),   # The directory you want to save the file in
    width = 14, # The width of the plot in inches
    height = 18) # The height of the plot in inches

figure

dev.off()


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sev_colors <- c("Mild" = "#3693a4","Severe" = "#f7464e")
singleR_colors <- c(
                    "CD4 T cells" = "#cecce2",
                    

                    "CD8 T cells" = "#422483",
                    "T regulatory cells" = "#2907b4",

                    
                    "gd T cells" = "#004c9d",

                    

                    "Natural killer cells" = "#338eb0",

                    
                    "MAIT cells" = "#d9dada",
                    

                    "B cells" = "#00963f",

                    "Plasmablasts" = "#d5e7dd",
                    


                    "Plasmacytoid dendritic cells" = "#ef7c00",
                    
                    "Myeloid dendritic cells" = "#e2a9cd",
                    
                    "Intermediate monocytes" = "#e6330f",
                    "Non classical monocytes" = "#ea5552",

                    "Classical monocytes" = "#f4a5a5",

                    

                    "Low-density neutrophils" = "#87cbbe",

                    "Platelet cells" = "#2a3937",
                    
                    "mixed"="#63CDE3"
)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

options(repr.plot.width=20, repr.plot.height=20)
c25 <- c(
  "#E69F00", "#56B4E9","#ed07ac" , "#F0E442", 
                            "#0072B2", "#D55E00", "#CC79A7", "#000000",
  "blue", "#E31A1C", # red
  "dodgerblue2",
  "#009E73", "green1", "yellow4"
  ,"#6707ed"
)

annotate <- function(bonn){
Idents(bonn)=factor(bonn$celltype, levels = identities)
return(bonn)
}
bonn =annotate (bonn)
stan =annotate (stan)
kor =annotate (kor)
berlin =annotate (berlin)
p1 <- DimPlot(object = bonn, reduction = 'umap', label = FALSE, cols= singleR_colors,raster=FALSE)+ NoLegend()+theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"),axis.text=element_blank(),        axis.ticks=element_blank())
p2 <- DimPlot(object = berlin, reduction = 'umap', label = FALSE,cols= singleR_colors) + NoLegend()+theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"),axis.text=element_blank(),        axis.ticks=element_blank())
p3 <- DimPlot(object = kor, reduction = 'umap', label = FALSE,cols= singleR_colors) + 
   theme(legend.text = element_text(size=8))+theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"),axis.text=element_blank(),        axis.ticks=element_blank())
p4 <- DimPlot(object = stan, reduction = 'umap', label = FALSE,cols= singleR_colors,raster=FALSE) + NoLegend()+theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"),axis.text=element_blank(),        axis.ticks=element_blank())


Idents(bonn)=factor(bonn$condition)
Idents(berlin)=factor(berlin$condition)
Idents(kor)=factor(kor$condition)
Idents(stan)=factor(stan$condition)
p5 <- DimPlot(object = bonn, reduction = 'umap', label = FALSE, cols= sev_colors,raster=FALSE)+ NoLegend()+theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"),axis.text=element_blank(),        axis.ticks=element_blank())
p6 <- DimPlot(object = berlin, reduction = 'umap', label = FALSE,cols= sev_colors) + NoLegend()+theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"),axis.text=element_blank(),        axis.ticks=element_blank())

p7 <- DimPlot(object = kor, reduction = 'umap', label = FALSE,cols= sev_colors) + 
   theme(legend.text = element_text(size=8))+theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"),axis.text=element_blank(),        axis.ticks=element_blank())
p8 <- DimPlot(object = stan, reduction = 'umap', label = FALSE,cols= sev_colors,raster=FALSE) + NoLegend()+theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"),axis.text=element_blank(),        axis.ticks=element_blank())


figure <- ggarrange(p5, p1, p6,p2,p8,p4,p7,p3,
                    font.label = list(size = 16, color = "black", face = "bold", family = NULL),
                    labels = c("Bonn","Bonn", "Berlin", "Berlin","Stanford","Stanford",'Korea','Korea'),vjust= 1.2,
                    ncol = 2, nrow = 4)
figure

pdf(file = paste0("~/new_NUM_CODEX_PLUS/Amina/CellSubmission/new_runs/figures/umap2.pdf"),   # The directory you want to save the file in
    width = 14, # The width of the plot in inches
    height = 17) # The height of the plot in inches

figure

dev.off()


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bonn_orig <- readRDS("~/new_NUM_CODEX_PLUS/Covid19_data/cohort2.annote.rds")
berlin_orig <- readRDS("~/new_NUM_CODEX_PLUS/Covid19_data/cohort1.annote.rds")
Stan_orig <- readRDS("~/new_NUM_CODEX_PLUS/Amina/data_/blish_awilk_covid_seurat.rds")
identities_org= c("CD14+ Monocytes" ,"CD16+ Monocytes" ,"pDCs" , "mDCs" ,  "B"  ,  "Plasmablasts"   , "CD8+ T"   , "CD4+ T" , "Prol. T" ,"NK" ,  "Neutrophils","Immature Neutrophils",     "mix/undefined"  ,  "Megakaryocytes"    )
Idents(bonn_orig)=factor(bonn_orig$celltypeL0, levels = identities_org)
Idents(berlin_orig)=factor(berlin_orig$celltypeL0, levels = identities_org)

#DotPlot(bonn_orig , features = feat, cols = 'RdBu')+RotatedAxis()+ggtitle('Bonn Original Celltypes')+ labs(y="Split by Celltype", x = "Marker Genes")

#DotPlot(berlin_orig, features = feat, cols = 'RdBu')+RotatedAxis()+ggtitle('Berlin Original Celltypes')+ labs(y="Split by Celltype", x = "Marker Genes")
identities_org= c("CD14 Mono","CD16 Mono" ,"pDC" ,"cDC1","cDC2", "B naive","B memory","B intermediate","Plasmablast" ,"CD8 Naive" ,"CD8 TEM","CD8 TCM","CD4 Naive","CD4 TCM","CD4 TEM","CD4 CTL" ,"Treg","MAIT","gdT","dnT","NK","NK_CD56","Neutrophil","Developing neutrophil","HSPC","Proliferating","Eryth","ASDC","Platelet" )
Stan_orig$celltypeL0= Idents(Stan_orig)
Stan_orig = subset(Stan_orig, subset = celltypeL0 %in% c("CD14 Mono","CD16 Mono" ,"pDC" ,"cDC1","cDC2", "B naive","B memory","B intermediate","Plasmablast" ,"CD8 Naive" ,"CD8 TEM","CD8 TCM","CD4 Naive","CD4 TCM","CD4 TEM","CD4 CTL" ,"Treg","MAIT","gdT","dnT","NK","NK_CD56","Neutrophil","Developing neutrophil","HSPC","Platelet"))

Idents(Stan_orig)=factor(Stan_orig$celltypeL0, levels = identities_org)

#DotPlot(Stan_orig, features = feat, cols = 'RdBu')+RotatedAxis()+ggtitle('Stanford Original Celltypes')+ labs(y="Split by Celltype", x = "Marker Genes")



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bonn =annotate (bonn)
stan =annotate (stan)
kor =annotate (kor)
berlin =annotate (berlin)
feat=c('CD14','LYZ','VCAN','FCGR3A','LST1','MS4A7','APOBEC3A','IRF4','ITM2C','HLA-DQA1','CD1C','CLEC10A','HLA-DPA1','HLA-DRB1','HLA-DRA','CD19','MS4A1','PAX5',
       'MZB1' ,'POU2AF1','SLAMF7','PRDM1','TNFRSF17','TRAC', 'CD3G','CD3E','CD3D', 'CCR7','CCL5','CD8A','CD8B',
       'IL7R', 'LTB', 'LDHB', 'TPT1',  'TMSB10', 'KLRB1', 'GZMK','TRGC2', 'KLRD1', 'GZMH','GZMA','CST7','FGFBP2','NKG7','IL2RB','GNLY','FCGR3B','MME', 'PPBP','PF4','NRGN')
#feat=c(	'CD14','LYZ','VCAN','KIT', 'TRDC', 'TTLL10', 'LINC01229', 'SOX4', 'KLRB1', 'TNFRSF18', 'TNFRSF4', 'IL1R1', 'HPGDS','SPINK2', 'PRSS57', 'CYTL1', 'EGFL7', 'GATA2', 'CD34', 'SMIM24', 'AVP', 'MYB', 'LAPTM4B','HBD', 'HBM', 'AHSP', 'ALAS2', 'CA1', 'SLC4A1', 'IFIT1B', 'TRIM58', 'SELENBP1', 'TMCC2','PPBP', 'PF4', 'NRGN', 'GNG11', 'CAVIN2', 'TUBB1', 'CLU', 'HIST1H2AC', 'RGS18', 'GP9')
p1=DotPlot(bonn , features = feat, cols = 'RdBu')+RotatedAxis()+theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"),axis.text.y=element_text(hjust=1, size=8),        axis.text.x=element_text(angle=90, hjust=1, size=8),axis.title=element_blank())+NoLegend()

p2=DotPlot(berlin, features = feat, cols = 'RdBu')+theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"),axis.text.y=element_text(hjust=1, size=8),        axis.text.x=element_text(angle=90, hjust=1, size=8),axis.title=element_blank())+NoLegend()

p4= DotPlot(kor, features = feat, cols = 'RdBu')+theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"),axis.text.y=element_text(hjust=1, size=8),        axis.text.x=element_text(angle=90, hjust=1, size=8),axis.title=element_blank())

p3= DotPlot(stan, features = feat, cols = 'RdBu')+theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"),axis.text.y=element_text(hjust=1, size=8),        axis.text.x=element_text(angle=90, hjust=1, size=8),axis.title=element_blank())+NoLegend()

p5=DotPlot(bonn_orig , features = feat, cols = 'RdBu')+theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"),axis.text.y=element_text(hjust=1, size=8),        axis.text.x=element_text(angle=90, hjust=1, size=8),axis.title=element_blank())+NoLegend()

p6=DotPlot(berlin_orig, features = feat, cols = 'RdBu')+theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"),axis.text.y=element_text(hjust=1, size=8),        axis.text.x=element_text(angle=90, hjust=1, size=8),axis.title=element_blank())+NoLegend()

p7= DotPlot(Stan_orig, features = feat, cols = 'RdBu')+theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"),axis.text.y=element_text(hjust=1, size=8),        axis.text.x=element_text(angle=90, hjust=1, size=8),axis.title=element_blank())+NoLegend()
leg <- get_legend(p4)
p4 <- p4 + theme(legend.position = "none")
figure <- ggarrange(p1,p5,p2,p6,p3,p7,p4,leg,
                    font.label = list(size = 16, color = "black", face = "bold", family = NULL),
                    labels = c("A.1","A.2", "B.1", "B.2","C.1","C.2",'D'),vjust= 1.2,
                    ncol = 2, nrow = 4)

pdf(file = paste0("~/new_NUM_CODEX_PLUS/Amina/CellSubmission/new_runs/figures/Kownmarkergenes2.pdf"),   # The directory you want to save the file in
    width = 16, # The width of the plot in inches
    height = 20) # The height of the plot in inches

figure

dev.off()



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# DotPlot(bonn,split.by = 'condition', features = feat, cols = 'RdBu')+RotatedAxis()+ggtitle('Bonn')+ labs(y="Split by Celltype and Condition", x = "Marker Genes")
# DotPlot(berlin,split.by = 'condition', features = feat, cols = 'RdBu')+RotatedAxis()+ggtitle('Berlin')+ labs(y="Split by Celltype and Condition", x = "Marker Genes")
# DotPlot(kor,split.by = 'condition', features = feat, cols = 'RdBu')+RotatedAxis()+ggtitle('Korean')+ labs(y="Split by Celltype and Condition", x = "Marker Genes")
# DotPlot(stan,split.by = 'condition', features = feat, cols = 'RdBu')+RotatedAxis()+ggtitle('Stanford')+ labs(y="Split by Celltype and Condition", x = "Marker Genes")




## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Idents(bonn) = bonn$monaco.labels
# Idents(berlin) = berlin$monaco.labels
# Idents(stan) = stan$monaco.labels
# Idents(kor) = kor$monaco.labels
# feat=c('CD14','LYZ','VCAN','FCGR3A','LST1','MS4A7','APOBEC3A','IRF8','ITM2C','CD1C','CLEC10A','HLA-DPA1','HLA-DRB1','HLA-DRA','CD19','MS4A1','PAX5','CR2','FCER2','CXCR5','TNFRSF13B','CD40','SLAMF7','PRDM1','CD38','CD27','TNFRSF17','TRAC', 'CD3G','CD3E','CD3D','CD4','CD8A','CD8B','CCL5','CCL4','GZMH','GZMA','CST7', 'NPM1', 'CCR7','KLRB1','FGFBP2','NKG7','NCAM1','IL2RB','GNLY','FCGR3B','MME','ITGAM','CCR3','CD69','PTGDR2','CD34','ALDH1A1')
# 
# DotPlot(bonn , features = feat, cols = 'RdBu')+RotatedAxis()+ggtitle('Bonn')+ labs(y="Split by Celltype", x = "Marker Genes")
# 
# DotPlot(berlin, features = feat, cols = 'RdBu')+RotatedAxis()+ggtitle('Berlin')+ labs(y="Split by Celltype", x = "Marker Genes")
# 
# DotPlot(kor, features = feat, cols = 'RdBu')+RotatedAxis()+ggtitle('Korean')+ labs(y="Split by Celltype", x = "Marker Genes")
# 
# DotPlot(stan, features = feat, cols = 'RdBu')+RotatedAxis()+ggtitle('Stanford')+ labs(y="Split by Celltype", x = "Marker Genes")



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ggplot(data = as.data.frame(table(bonn$monaco.labels.f)), aes(x=Var1, y=Freq)) + geom_bar(stat = "identity")+ coord_flip()
# ggplot(data = as.data.frame(table(berlin$monaco.labels.f)), aes(x=Var1, y=Freq)) + geom_bar(stat = "identity")+ coord_flip()
# ggplot(data = as.data.frame(table(stan$monaco.labels.f)), aes(x=Var1, y=Freq)) + geom_bar(stat = "identity")+ coord_flip()
# ggplot(data = as.data.frame(table(kor$monaco.labels.f)), aes(x=Var1, y=Freq)) + geom_bar(stat = "identity")+ coord_flip()



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

sessionInfo()

