info = read.csv("OutputRCC1.csv")
rownames(info) = info[,1]
info = info[,-1]

num = grep("Level_", colnames(info))

finalOut  = "."

for(i in 1:length(num)){
  j = num[i]
  if(colnames(info)[j] == "Level_1"){
    clusterAnnotation("OutputRCC1.csv",  "cellType_new", "Level_1", finalOut)
  } else {
    cols = num[1:i]
    info$Clusters_new = apply(info[ ,cols], 1, paste, collapse = "_")
    colnames(info)[ncol(info)] = paste0(colnames(info)[num[i]], "new")
    write.csv(info, file = "clusterAnno.csv")
    clusterAnnotation("clusterAnno.csv",  "cellType_new", paste0(colnames(info)[num[i]], "new"), finalOut)
  }
}
dev.off()
