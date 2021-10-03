infoFile = "ligerOut.csv"
data = readRDS("log2SNU.rds")


info = read.csv(infoFile)
rownames(info) = info[,1]
info = info[,-1]

identical(colnames(data), rownames(info))

data = data[, order(colnames(data))]
info = info[order(rownames(info)), ]

a = findMarkers(as.matrix(data), info$seurat_clusters, direction = "up", lfc = 1)
mat = NULL

for(i in 1:length(a)){
  print(i)
  clusterInfo = a[[i]]
  clusterInfo = subset(clusterInfo, clusterInfo$FDR < 0.01)
  lfc_mat = as.data.frame(clusterInfo[,4:ncol(clusterInfo)])
  lfc = apply(lfc_mat, 1, min)
  clusterInfo = cbind(clusterInfo, lfc)
  clusterInfo = subset(clusterInfo, clusterInfo$lfc > 1)
  clusterInfo = subset(clusterInfo, clusterInfo$lfc > 1)
  
  if(nrow(clusterInfo) == 0){
    mat = qpcR:::cbind.na(mat, rep("0", 5))
  } else {
    mat = qpcR:::cbind.na(mat, rownames(clusterInfo))
  }
}
