clusters = unique(info$Clusters)

cells = unique(info$cellType)
samples = unique(info$sampleID)
 
cols = append(cells, samples)
cols = append(cells, "Total")

mat = matrix(0, ncol = length(cols),  nrow = length(clusters))
rownames(mat) = clusters
colnames(mat) = cols

for(i in 1:length(clusters)){
  sub = subset(info, info$Clusters == clusters[i])
  mat[as.character(clusters[i]), "Total"] = nrow(sub)
  a = as.data.frame(table(sub$cellType))
  for(j in 1:nrow(a)){
    mat[as.character(clusters[i]), as.character(a[j,1])] = a[j,2]
  }
  
  #a = as.data.frame(table(sub$sampleID))
  #for(j in 1:nrow(a)){
    #mat[as.character(clusters[i]), as.character(a[j,1])] = a[j,2]
  #}
}


pdf("cellHighlight1.pdf", width = 20, height = 20)

cellTypes = unique(data@meta.data$cellType)

for(i in 1:length(cellTypes)){
  num = which(data@meta.data$cellType == as.character(cellTypes[i]))
  cellNames = rownames(data@meta.data)[num]
  p1 = DimPlot(object = data, cells.highlight = cellNames, cols.highlight = "red", 
               cols = "gray", order = F, reduction = "numap", label = T, label.size = 1, 
               pt.size = 0.01, sizes.highlight = 0.01) + plot_annotation(title = cellTypes[i])
  print(p1)
}

dev.off()
