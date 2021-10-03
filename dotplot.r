data = readRDS("../log2SNU.rds")

info = read.csv("OutputRCC1.csv")
rownames(info) = info[,1]
info = info[,-1]

markers = read.csv("finalClusters.csv")

obj = readRDS("Level1.rds")
identical(rownames(info), rownames(obj@meta.data))

obj@meta.data = info

library(Seurat)
library(ggplot2)
library(ggpubr)

levels = unique(markers[,2])
cols = grep("Level_", colnames(obj@meta.data))
cc = grep("cellType_l1",colnames(obj@meta.data))

uPlots = list()
genesF = NULL
'%ni%' <- Negate('%in%')

for(i in 1:length(levels)){
  imarkers = subset(markers, markers[,2] == i)
  if(i > 1){
    imarkers = subset(imarkers, imarkers[,1] %ni% genesF)
  }
  tmarkers = imarkers %>% group_by(clusters) %>% dplyr::slice(1:3)
  tmarkers = as.data.frame(tmarkers)
  genes = tmarkers[,1]
  genes = unique(genes)
  genesF = append(genesF, genes)
  genesF = unique(genesF)
  if(length(genes) < 1){  next }
  if(i == 1){
   obj@active.ident = factor(obj@meta.data$Level_1)
   names(obj@active.ident) = rownames(obj@meta.data)
   p1 = DotPlot(object = obj,  features = genesF) + ggtitle("Level_1") +  RotatedAxis()
  } else {
    cols1 = cols[1:i]
    #cols1 = append(cc, cols1)
    meta = obj@meta.data
    meta$newClus = apply(meta[ ,cols1], 1, paste0, collapse = "_")
    meta$newClus = gsub(" ", "", meta$newClus)
    obj@active.ident = factor(meta$newClus)
    names(obj@active.ident) = rownames(meta)
    p1 = DotPlot(object = obj,  features = genesF) + ggtitle(paste0("Level_", i)) +  RotatedAxis()
  }
  uPlots[[i]] = p1
}

pdf("dotplot_Levels1.pdf", width = 100, height = round(length(uPlots)/2, 0) * 30)
ggarrange(plotlist = uPlots, ncol = 2, nrow = round(length(uPlots)/2, 0))
dev.off()

