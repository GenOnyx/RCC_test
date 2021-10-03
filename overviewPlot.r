
######################################
data = readRDS("finalSeuratObj.rds")

info = read.csv("OutputRCC1.csv")
rownames(info) = info[,1]
info = info[,-1]

markers = read.csv("LevelMarkers/combineMarkers.csv")

cpm = readRDS("../log2SNU.rds")

identical(rownames(info), colnames(cpm))
identical(rownames(info), rownames(data@meta.data))

data@meta.data = info
identical(rownames(data@meta.data), rownames(data@reductions$numap@cell.embeddings))

info = cbind.data.frame(info, data@reductions$numap@cell.embeddings)
clus = unique(info$Clusters)
sample_percent = as.data.frame(table(info$sample_location))
rownames(sample_percent) = sample_percent[,1]
colnames(sample_percent) = c("SL", "freq")

for(i in 1:length(clus)){
  num = which(info$Clusters == i)
  #cells = rep("grey", nrow(info))
  #cells[num] = "red"
  
  i2 = cbind.data.frame(info, cells)
  i1 = subset(info, info$Clusters == i)
  sl = as.data.frame(table(i1$sample_location))
  sl[,3] = (sl[,2]/sum(sl[,2]))*100
  sl[,4] = as.character(sl[,1])
  sl[,4] = sl[,4] %>%
    gsub("Core_SNU27_Core", "lightblue", .) %>%
    gsub("Core_SNU33_Core", "pink", .) %>% 
    gsub("Core_SNU33b2_Core", "hotpink", .) %>%
    gsub("Core_SNU40_Core", "springgreen", .) %>%
    gsub("Core_SNU43_Core", "yellow", .) %>%
    gsub("Accutase_UK", "purple", .) %>%
    gsub("GBM_v7_UK", "mediumorchid4", .) %>%
    gsub("Peri_SNU27_Peri", "navy", .) %>%
    gsub("Peri_SNU33_Peri", "hotpink3", .) %>%
    gsub("Peri_SNU33b2_Peri", "hotpink4", .) %>%
    gsub("Peri_SNU40_Peri", "olivedrab", .) %>%
    gsub("Peri_SNU43_Peri", "goldenrod1", .) %>%
    gsub("SNU38_Control", "black", .)
  
  colnames(sl) = c("location", "val", "percent", "color")
  rownames(sl) = sl[,1]
  sl = merge(sl, sample_percent, all.x = T, by = 0)
  sl = sl[,-1]
  sl$spercent = (sl$val/sl$freq) * 100
  sl$label = paste0(sl$val, "(", round(sl$spercent, 1), "%)")
  
  
  fclus = unique(i1$l6)
  fclus = paste0("K", fclus)
  features = as.vector(markers[, fclus])
  features = gsub(0, NA, features)
  features = as.vector(na.omit(unique(features)))
  
  if(length(features) == 0){ next }
  
  pdf(paste0("overview_", fclus, ".pdf"), width = 40, height = 20)
  
  m = markers[,fclus]
  m = unique(m)
  m = as.vector(na.omit(m))
  if(length(m) > 5){
    m = m[1:5]
  }
  #p1 = plot_density(data, reduction = "numap", features = m, joint = T, combine = T)
  #sl$clus = rep(fclus, nrow(sl))
  
  par(mfrow = c(1,1))
  if(length(m) > 1){
    mat = subset(cpm, rownames(cpm) %in% m)
    mat = t(mat)
    meta = cbind.data.frame(info, mat)
    rgb_plot(meta, m)
  } else{
    mat = subset(cpm, rownames(cpm) %in% m)
    mat = t(mat)
    meta = cbind.data.frame(info, mat)
    fp = ggplot(meta, aes_string(x = "UMAP_1", y = "UMAP_2", col = m[1])) + 
      geom_point(size = 0.01)  + theme_classic() +
      theme(plot.title = element_text(size = 30)) + 
      scale_color_gradientn(colors = c("lightgrey", "red"))
    h = layer_data(fp)
    plot(h$x, h$y, pch = 20, cex = 0.1, 
         xlab = fclus, ylab = "", main = m[1], 
         col = h$colour)
  }
  
  convexHull(i2, fclus, paste0("K",clus[i]))
  
  #par(mfrow = c(1,2))
  pie(sl$percent, sl$label, col = sl$color, cex = 5)
  legend("topright", legend = sl$location, col = sl$color, 
         pch = 20, cex = 1)
 
  plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
       xaxt='n', yaxt='n', xlab='', ylab='', main = "Markers")
  pos1 = 0.1
  for(pos in length(features):1){
    text(1, pos1, features[pos], pos = 1, cex = 5)
    pos1 = pos1 + 0.1  
  }
  
  if(length(m) == 1){
    mark1 = plot_density(data, reduction = "numap", features = m) + NoLegend()
    mark2 = 
    print(ggarrange(p1[[length(m) + 1]],p2,ncol = 2, nrow = 1))
    
  }
  
  #if(length(m) > 1){
    #print(ggarrange(p1[[length(m) + 1]],p2,ncol = 2, nrow = 1))
    #dev.off()
  #} else {
    #print(ggarrange(p1,p2,ncol = 2, nrow = 1))
    #dev.off()
  #}
}

