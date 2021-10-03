LevelHeatmap = function(expr, info, genes, attr1){
  library(ComplexHeatmap)
  library(circlize)
  library(matrixStats)
  
  data = readRDS(expr)
  data = data[,  order(colnames(data))]
  
  info = read.csv(info)
  rownames(info) = info[,1]
  info = info[,-1]
  info = info[order(rownames(info)), ]
  identical(colnames(data),  rownames(info))
  
  markers = read.csv(genes)
  markers = unique(markers[,1])
  data = subset(data, rownames(data) %in% markers)
  
  cmat = cbind.data.frame(info$concatenate, t(data))
  
  ctypes = unique(info[, attr1])

  pdf("Heatmap_markers.pdf", width = 40, height = 50)
  for(i in 1:length(ctypes)){
    i1 = subset(info, info[, attr1] == ctypes[i])
    num = which(colnames(data) %in% rownames(i1))
    mat = data[,num]
    a = apply(mat, 1, sum)
    a = subset(a, a > 0)
    mat = subset(mat, rownames(mat) %in% names(a))
    identical(colnames(mat), rownames(i1))
    mat2 = cbind.data.frame(i1[, attr2], t(mat))
    mat2 = aggregate(mat2[, 2:ncol(mat2)], list(mat2[,1]), mean)
    rownames(mat2) = mat2[,1]
    mat2 = mat2[,-1]
    mat2 = t(mat2)
  
    mat1 = (mat2 - rowMeans(mat2))/(rowSds(as.matrix(mat2)))[row(mat2)]
    print(Heatmap(mat1, 
            col = colorRamp2(c(-2,0,2), c("orangered", "white",  "purple")), 
            column_names_gp = gpar(fontsize = 3),
            row_names_gp = gpar(fontsize = 0),  column_title = ctypes[i], 
            column_title_gp = gpar(fontsize = 50)))
  }
  dev.off()
}
