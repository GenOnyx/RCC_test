
setwd("LeidenRCC/")
data = readRDS("Level1.rds")
umap = read.csv("UMAPcoords.csv") 
rownames(umap) = umap[,1]
umap = umap[,-1]
colnames(umap)= c("UMAP_1", "UMAP_2")
info = read.csv("OutputRCC1.csv")
rownames(info) = info[,1]
info = info[,-1]


identical(rownames(data@reductions$umap@cell.embeddings), rownames(umap))
data@reductions$numap = data@reductions$umap
data@reductions$numap@cell.embeddings = as.matrix(umap)

identical(rownames(data@meta.data), rownames(info))
data@meta.data = info

data@active.ident = factor(data@meta.data$Level_1)
names(data@active.ident) = rownames(data@meta.data)

genes = read.csv("finalClusters.csv")
num = which(genes[,2] == 2)
genes1 = genes[num,]

features = genes1 %>% group_by(clusters) %>% slice(1:10)
features = as.data.frame(features)
features = features[,1]

p1 = DimPlot(data,  label = T, group.by = "Clusters", reduction = "umap")
p2 = DimPlot(data,  label = T, group.by = "Clusters", reduction = "numap")
ggarrange(p1,p2, ncol = 2)
