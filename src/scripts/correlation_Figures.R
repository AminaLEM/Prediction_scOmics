library(ggplot2)
library(tidyverse)
library(Seurat)
library(ggpubr)

set.seed(0)


path_to_outData= snakemake@params[[1]]
path_to_outFig= snakemake@params[[2]]


plot_correlation = function(bonn,sam_name){
bonn_cells = as.data.frame(cbind(colnames(bonn),bonn$sampleID, bonn$celltype, bonn$condition))
colnames(bonn_cells)= c('cells', 'sampleID', 'celltypes','condition')
coun= count(bonn_cells %>% group_by(sampleID, condition))
bonn_cells = bonn_cells[bonn_cells$sampleID %in% coun[coun['n']>1000,'sampleID']$sampleID,]
print(length(unique(bonn_cells$sampleID)))
print(count(bonn_cells %>% group_by(sampleID, condition)))
for (ct in unique(bonn$celltype)){
  #ct= 'Classical monocytes'
  bonn_cell = bonn_cells[bonn_cells$celltypes==ct,]
  min_cells = min(count(bonn_cell %>% group_by(sampleID))[,'n'])    
  print(min_cells)
  bonn_cell_sampled = bonn_cell %>% group_by(sampleID)%>% slice_sample(n=min_cells)
  if (ct == unique(bonn$celltype)[1]){ sampled_cells =bonn_cell_sampled
  }else{
  sampled_cells = rbind(sampled_cells,bonn_cell_sampled)
  }
}
seu.filtered =  bonn[,colnames(bonn) %in% sampled_cells$cells]
seu.filtered$identif <-
  paste(seu.filtered@meta.data$sampleID,
        seu.filtered@meta.data$condition,
        seu.filtered@meta.data$who_score,
        seu.filtered@meta.data$batch,        
        sep = "__")

ch_avr_nor <- AggregateExpression(seu.filtered, 
                                  group.by = c("identif"),
                                  assays = 'RNA',
                                  slot = "data",
                                  return.seurat = TRUE)

ch_avr_nor@meta.data[c('sampleID', 'condition', 'who_score','batch')] = t(data.frame(strsplit(colnames(ch_avr_nor), "__")))

bonn$identif <-
  paste(bonn@meta.data$sampleID,
        bonn@meta.data$condition,
        bonn@meta.data$who_score,
        bonn@meta.data$batch,        
        sep = "__")
subset_mono = subset(bonn, subset= sampleID %in% coun[coun['n']>1000,'sampleID']$sampleID)
bonn_avr_nor <- AggregateExpression(subset_mono, 
                                  group.by = c("identif"),
                                  assays = 'RNA',
                                  slot = "data",
                                  return.seurat = TRUE)
bonn_avr_nor@meta.data[c('sampleID', 'condition', 'who_score','batch')] = t(data.frame(strsplit(colnames(bonn_avr_nor), "__")))

print(cor(ch_avr_nor@assays$RNA@data['RETN',], bonn_avr_nor@assays$RNA@data['RETN',], method =meth))
print(cor(ch_avr_nor@assays$RNA@data['S100P',], bonn_avr_nor@assays$RNA@data['S100P',], method = meth))
print(cor(ch_avr_nor@assays$RNA@data['ANXA3',], bonn_avr_nor@assays$RNA@data['ANXA3',], method = meth))
my_data = as.data.frame(cbind(ch_avr_nor$condition,ch_avr_nor@assays$RNA@data['RETN',],bonn_avr_nor@assays$RNA@data['RETN',]))
print(my_data)
colnames(my_data)= c('condition',"subsampled","all_cells")

my_data$subsampled <- as.numeric(as.character(my_data$subsampled) )
my_data$all_cells <- as.numeric(as.character(my_data$all_cells) )

p1=ggscatter(my_data, x = "subsampled", y = "all_cells",color = 'condition',
          add = "reg.line", conf.int = TRUE, 
          xlab = "", ylab = "normalized expression by all cells",palette = c("#3693a4", "#f7464e"))
  p1=p1+  stat_cor(size = 7, method = meth)+
  theme(legend.position = "none",plot.margin = margin(0.8,0.5,0.5,0.5, "cm"),text = element_text(size = 15))

my_data = as.data.frame(cbind(ch_avr_nor$condition,ch_avr_nor@assays$RNA@data['S100P',],bonn_avr_nor@assays$RNA@data['S100P',]))
print(my_data)
colnames(my_data)= c('condition',"subsampled","all_cells")
my_data$subsampled <- as.numeric(as.character(my_data$subsampled) )
my_data$all_cells <- as.numeric(as.character(my_data$all_cells) )

p2=ggscatter(my_data, x = "subsampled", y = "all_cells", color = 'condition',,
             add = "reg.line", conf.int = TRUE, 
             xlab = "normalized expression by subsampled cells", ylab = "",palette = c("#3693a4", "#f7464e"))
p2=p2+  stat_cor(size = 7, method = meth)+
  theme(legend.position = "none",plot.margin = margin(0.8,0.5,0.5,0.5, "cm"),text = element_text(size = 15))

my_data = as.data.frame(cbind(ch_avr_nor$condition,ch_avr_nor@assays$RNA@data['ANXA3',],bonn_avr_nor@assays$RNA@data['ANXA3',]))
print(my_data)
colnames(my_data)= c('condition',"subsampled","all_cells")

my_data$subsampled <- as.numeric(as.character(my_data$subsampled) )
my_data$all_cells <- as.numeric(as.character(my_data$all_cells) )


p3=ggscatter(my_data, x = "subsampled", y = "all_cells", color = 'condition',
             add = "reg.line", conf.int = TRUE, 
             xlab = "", ylab = "",palette = c("#3693a4", "#f7464e"))
  p3=p3+  stat_cor(size = 7, method = meth)+
  theme(legend.position = "none",plot.margin = margin(0.8,0.5,0.5,0.5, "cm"),text = element_text(size = 15))
figure <- ggarrange(p1,p2,p3,
                    font.label = list(size = 20, color = "black", face = "bold", family = NULL),
                    labels = c("RETN", "S100P","ANXA3"),vjust= 1.2,
                    nrow = 1)
pdf(file = paste0(paste0(paste0(paste0(path_to_outFig,"/"),meth),sam_name),"correlation.pdf"),   # The directory you want to save the file in
    width = 14, # The width of the plot in inches
    height = 5) # The height of the plot in inches

print(figure)

dev.off()
}

bonn= readRDS(paste0(path_to_outData,'/bonn_anno.rds'))
berlin= readRDS(paste0(path_to_outData,'/berlin_anno.rds'))
kor= readRDS(paste0(path_to_outData,'/kor_anno.rds'))
stan= readRDS(paste0(path_to_outData,'/stan_anno.rds'))

training_set= merge(bonn, y = berlin, add.cell.ids = c("bonn", "berlin"), project = "training_set",merge.data = TRUE)
testing_set= merge(stan, y = kor, add.cell.ids = c("stan", "kor"), project = "testing_set",merge.data = TRUE)

meth="pearson"
plot_correlation(training_set,'train')
plot_correlation(testing_set,'test')

meth="spearman"
plot_correlation(training_set,'train')
plot_correlation(testing_set,'test')


plot_correlation(berlin,'berlin')
plot_correlation(stan,'stan')

