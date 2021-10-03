UMAPplot_func = function(LevMarkers, m1){
library(Seurat)
library(ggplot2)
library(ggpubr)
library(Nebulosa)

source("/media/user/disk3/R/convexHull_func.r")
source("/media/user/disk3/R/rgb_plot.r")
source("/media/user/disk3/R/piechart_func.r")

data = read.csv("UMAPcoords.csv")
rownames(data) = data[,1]
data = data[,-1]
colnames(data) = gsub("parent", "", colnames(data))


info = read.csv("OutputRCC1.csv")
rownames(info) = info[,1]
info = info[,-1]

cpm = readRDS("log2SNU.rds")

identical(rownames(data), rownames(info))

info = cbind.data.frame(info, data)

Tumap = read.csv("Tumap.csv")
rownames(Tumap) = Tumap[,1]
Tumap = Tumap[,-1]
print("file done")

sample_percent = as.data.frame(table(info$sample_location))
rownames(sample_percent) = sample_percent[,1]
colnames(sample_percent) = c("SL", "freq")

for(i in 1:length(LevMarkers)){
  fclus = LevMarkers[i]
  print(fclus)
  lev = paste0("Level_", length(strsplit(fclus, "_")[[1]]))
  
  nclus = fclus %>%
    gsub("K", "", .) %>%
    gsub("_", "", .)
  nLev = gsub("Level_", "", lev)
  
  if(lev == "Level_1"){
    i1 = info
    fclus1 = fclus
  } else {
    fclus1 = gsub("_0", "", fclus)
    lev1 = paste0("l", length(strsplit(fclus1, "_")[[1]]) - 1)
    if(lev1 == "l1"){
      lev1 = "Level_1"
    }
    split = strsplit(fclus1, "_")[[1]]
    split = split[-length(split)]
    nclus1 = paste0(split, collapse="_")
    nclus1 = gsub("K", "", nclus1)
    i1 = subset(info, info[, lev1]  == nclus1)
  }
  pdf(paste0(lev, "_", fclus, "_", paste(m1,collapse="_"), ".pdf"), width = 40, height = 40)

    mat = subset(cpm, rownames(cpm) %in% m1)
    meta1 = cbind.data.frame(info, t(mat))
    mat = t(mat)
    mat = subset(mat, rownames(mat) %in% rownames(i1))
    meta = cbind.data.frame(i1, mat)
    colnames(meta1) = gsub("-", "_", colnames(meta1))
    colnames(meta) = gsub("-", "_", colnames(meta))
    m1 = gsub("-", "_", m1)
    rgb_plot(meta1, m1, lev, fclus)
    convexHull(info, fclus1, lev, 0)
    rgb_plot(meta, m1, lev, fclus)
    convexHull(i1, fclus1, lev, 0)
  
    m1 = gsub("_", "-", m1)
    m2 = gsub("_", "-", m2)
  
    mat = subset(cpm, rownames(cpm) %in% m1)
    mat = t(mat)
    mat = subset(mat, rownames(mat) %in% rownames(i1))
    meta = cbind.data.frame(i1, mat)
    colnames(meta) = gsub("-", "_", colnames(meta))
    m1 = gsub("-", "_", m1)
      
    nLev = gsub("Level_", "", lev)
    nLev = as.numeric(nLev)
  
    num = grep("UMAP_", colnames(meta))
    print("col replaced")
    colnames(meta)[num] = c("parentUMAP_1", "parentUMAP_2")
    print("pre merge")
    meta = merge(meta, Tumap, by = 0, all.x = T)
    print("post merge")
    colnames(meta) = gsub(paste0("Lev", nLev, "_UMAP1"), "UMAP_1", colnames(meta))
    colnames(meta) = gsub(paste0("Lev", nLev, "_UMAP2"), "UMAP_2", colnames(meta))
      
    rownames(meta) = meta[,1]
    meta = meta[,-1]
      
    rgb_plot(meta, m1, lev, fclus)
    convexHull(meta, fclus1, nLev, 1)
    
    if(lev == "Level_1"){
      i2 = subset(info, info[,lev] == gsub("K", "", fclus))
    } else {
      i2 = subset(info, info[,gsub("Level_", "l", lev)] == gsub("K", "", fclus))
    }
    pval = piechart_func(i2, sample_percent)
  }
  
  dev.off()
}
