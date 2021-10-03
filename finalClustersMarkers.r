library(scran)
library(data.table)

DEG = function(fileName, infoFile){
  data = as.data.frame(fread(fileName))
  rownames(data) = data[,1]
  data = data[,-1]
  #data = t(data)
  
  info = read.csv(infoFile)
  rownames(info) = info[,1]
  info = info[,-1]
  genes = read.csv("genesUsed.csv")
  genes = genes[,1]
  
  mat = subset(data, rownames(data) %in% genes)
  mat = t(mat)
  ord_mat = mat[order(rownames(mat)), ]
  ord_info = info[order(rownames(info)),]
  
  identical(rownames(ord_mat), rownames(ord_info))
  ord_mat = t(ord_mat)
  Levels = which(grepl("Level_*", colnames(ord_info)))
  previousLevel = 0
  
  deg_info = ord_info
  
  level_markers = NULL
  
  for(li in 1:length(Levels)){
    colNums = Levels[li]
    print(colnames(ord_info)[colNums])
    if(li == 1){
      a = findMarkers(as.matrix(ord_mat), ord_info[,colNums], direction = "up", lfc = 1)
      mat = NULL
      clus = NULL
      for(i in 1:length(a)){
        print(i)
        clusterInfo = a[[i]]
        clusterInfo = subset(clusterInfo, clusterInfo$FDR < 0.01)
        lfc_mat = as.data.frame(clusterInfo[,4:ncol(clusterInfo)])
        lfc = apply(lfc_mat, 1, min)
        clusterInfo = cbind(clusterInfo, lfc)
        clusterInfo = subset(clusterInfo, clusterInfo$lfc > 1)
        print(dim(clusterInfo))
        if(nrow(clusterInfo) == 0){
          mat = append(mat, rep("0", 10))
        } else {
          mat = append(mat, rownames(clusterInfo)[1:10])
        }
        clus = append(clus, rep(names(a)[i], 10))
      }
      lev = rep(colnames(ord_info)[Levels[li]], length(mat))
      final = cbind.data.frame(mat,lev)
      final = cbind.data.frame(final, clus)
      write.csv(final, file = paste0("MH_", colnames(ord_info)[colNums], ".csv"))
      prev_li = li
    } else{
      clusters = unique(ord_info[, Levels[prev_li]])
      for(m in 1:length(clusters)){
        subInfo = subset(ord_info, ord_info[,colNums] == clusters[m])
        if(unique(ord_info[,colNums]) == 0){ next }
        mData = t(ord_mat)
        subData = subset(mData, rownames(mData) %in% rownames(subInfo))
        subData = t(subData)
        a = findMarkers(as.matrix(subData), subInfo[,colNums], direction = "up", lfc = 1)
        mat = NULL
        clus = NULL
        for(i in 1:length(a)){
          print(i)
          clusterInfo = a[[i]]
          clusterInfo = subset(clusterInfo, clusterInfo$FDR < 0.01)
          lfc_mat = as.data.frame(clusterInfo[,4:ncol(clusterInfo)])
          lfc = apply(lfc_mat, 1, min)
          clusterInfo = cbind(clusterInfo, lfc)
          clusterInfo = subset(clusterInfo, clusterInfo$lfc > 1)
          print(dim(clusterInfo))
          if(nrow(clusterInfo) == 0){
            mat = append(mat, rep("0", 10))
          } else {
            mat = append(mat, rownames(clusterInfo)[1:10])
          }
          mClus = apply(subInfo[,Levels[1]:colNums], 1, paste, collapse = "_")
          clus = append(clus, rep(names(a)[i], 10))
        }
        
      }
    }	
    
  }
}