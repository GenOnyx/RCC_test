library(Seurat)
library(ggplot2)
library(ggpubr)
library(Nebulosa)

source("/media/user/disk3/R/convexHull_func.r")
source("/media/user/disk3/R/piechart_func.r")
source("/media/user/disk3/R/rgb_plot.r")
source("/media/user/disk3/R/marker_print.r")

data = readRDS("finalSeuratObj.rds")

info = read.csv("OutputRCC1.csv")
rownames(info) = info[,1]
info = info[,-1]

LevMarkers = read.csv("LevelMarkers/combineMarkers.csv")
finMarkers = read.csv("finalClusters.csv")

cpm = readRDS("../log2SNU.rds")

identical(rownames(data@meta.data), rownames(info))

identical(rownames(data@reductions$numap@cell.embeddings), rownames(info))


info = cbind.data.frame(info, data@reductions$numap@cell.embeddings)
sample_percent = as.data.frame(table(info$sample_location))
rownames(sample_percent) = sample_percent[,1]
colnames(sample_percent) = c("SL", "freq")


Tumap = read.csv("Tumap.csv")
rownames(Tumap) = Tumap[,1]
Tumap = Tumap[,-1]
print("file done")



pmat = matrix(NA, nrow = length(LevMarkers), ncol = 2)

for(i in 1:ncol(LevMarkers)){
  fclus = colnames(LevMarkers)[i]
  print(fclus)
  lev = LevMarkers[1,i]
  m1 = LevMarkers[c(2:nrow(LevMarkers)), i]
  m1 = unique(m1)
  m1 = as.vector(na.omit(m1))
  if(length(m1) == 0){
    m1 = "0"
  }
  
  nclus = fclus %>%
    gsub("K", "", .) %>%
    gsub("_", "", .)
  nLev = gsub("Level_", "", lev)
  m2 = subset(finMarkers, (finMarkers[,2]== nLev) & (finMarkers[,3] == nclus))
  if(nrow(m2) == 0){
    m2 = "0"
  } else {
    m2 = m2[,1]
  }
  
  if((length(m1) == 1) & (length(m2) == 1)){
    if((m1 == "0") & (m2 == "0")){
      next
    }
  }
  
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
  pdf(paste0("overview1/", lev, "_", fclus, ".pdf"), width = 40, height = 40)
  if((length(m1) == 1) & (m1 == "0")){
    par(bg="gray40")
    plot(i1$UMAP_1, i1$UMAP_2, col = "gray40", main = "Final markers not present")
  } else {
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
  }
  if((length(m2) == 1) & (m2 == "0")){
    par(bg="gray40")
    plot(i1$UMAP_1, i1$UMAP_2, col = "gray40", main = "Level wise markers not present")
  } else {
    mat = subset(cpm, rownames(cpm) %in% m2)
    meta1 = cbind.data.frame(info, t(mat))
    mat = t(mat)
    mat = subset(mat, rownames(mat) %in% rownames(i1))
    meta = cbind.data.frame(i1, mat)
    colnames(meta1) = gsub("-", "_", colnames(meta1))
    colnames(meta) = gsub("-", "_", colnames(meta))
    m2 = gsub("-", "_", m2)
    
    rgb_plot(meta1, m2, lev, fclus)
    convexHull(info, fclus1, lev, 0)
    rgb_plot(meta, m2, lev, fclus)
    convexHull(i1, fclus1, lev, 0)
  }
  
  m1 = gsub("_", "-", m1)
  m2 = gsub("_", "-", m2)
  
  if((length(m1) > 0) & (length(m2) > 0)){
    if(m1 == "0"){
      print("useless")
    } else if(m2 == "0"){
      print("useless")
    } else {
      mat = subset(cpm, rownames(cpm) %in% m1)
      mat = t(mat)
      mat = subset(mat, rownames(mat) %in% rownames(i1))
      meta = cbind.data.frame(i1, mat)
      colnames(meta) = gsub("-", "_", colnames(meta))
      m1 = gsub("-", "_", m1)
      
      nLev = gsub("Level_", "", lev)
      nLev = as.numeric(nLev)
  
      #meta1 = meta
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
      
      par(mfrow = c(2,1))
      rgb_plot(meta, m1, lev, fclus)
      convexHull(meta, fclus1, nLev, 1)
     
      mat = subset(cpm, rownames(cpm) %in% m2)
      mat = t(mat)
      mat = subset(mat, rownames(mat) %in% rownames(i1))
      meta = cbind.data.frame(i1, mat)
      colnames(meta) = gsub("-", "_", colnames(meta))
      m2 = gsub("-", "_", m2)
     
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
      
      
      rgb_plot(meta, m2, lev, fclus)
      convexHull(meta, fclus1, nLev, 1)
       
    }
  }
  
  if(lev == "Level_1"){
    i2 = subset(info, info[,lev] == gsub("K", "", fclus))
  } else {
    i2 = subset(info, info[,gsub("Level_", "l", lev)] == gsub("K", "", fclus))
  }
  pval = piechart_func(i2, sample_percent)
  marker_print(m1,m2)
  dev.off()
  
  pmat[i, 1] = fclus
  pmat[i, 2] = as.numeric(pval)
}

 write.csv(pmat, file = "chisq.csv")
