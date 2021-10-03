DEG_clus = function(exprFile, infoFile, clus1, clus2){
  library(data.table)
  
  data = readRDS(exprFile)
  
  info = read.csv(infoFile)
  rownames(info) = info[,1]
  info = info[,-1]
  
  identical(colnames(data), rownames(info))
  
  lev1 = strsplit(clus1, "_")[[1]]
  lev1 = length(lev1)
  lev1 = paste0("l", lev1)
  if(lev1 == "l1"){
    lev1 = "Level_1"
  }
  
  lev2 = strsplit(clus2, "_")[[1]]
  lev2 = length(lev2)
  lev2 = paste0("l", lev2)
  if(lev2 == "l1"){
    lev2 = "Level_1"
  }
  
  i1 = subset(info, info[, lev1] == gsub("K", "", clus1))
  if(nrow(i1) == 0){
    return("Clus1 not present")
  }
  
  i2 = subset(info, info[, lev2] == gsub("K", "", clus2))
  if(nrow(i2) == 0){
    return("Clus2 not present")
  }
  
  num = which(colnames(data) %in% rownames(i1))
  mat1 = data[,num]
  
  num = which(colnames(data) %in% rownames(i2))
  mat2 = data[,num]
  
  mat = matrix(NA, ncol = 10, nrow = nrow(mat1))
  
  for(i in 1:nrow(mat)){
    mat[i,1] = rownames(mat1)[i]
    a = as.vector(mat1[i,])
    b = as.vector(mat2[i,])
    
    mat[i,2] = length(which(a > 0))
    mat[i,3] = length(which(b > 0))
    
    mat[i,4] = (as.numeric(mat[i,2])/ncol(mat1)) * 100
    mat[i,5] = (as.numeric(mat[i,3])/ncol(mat2)) * 100
    
    if((mean(a) == 0) & (mean(b) == 0)){next}
    
    mat[i,6] = mean(a)
    mat[i,7] = mean(b)
    
    res = t.test(a, b, paired = F, alternative = "two.sided")
    mat[i,8] = mean(a) - mean(b)
    mat[i,9] = as.numeric(res$p.value)
  }
  colnames(mat) = c("genes", "posCellsA", "posCellsB", 
                    "percentA", "percentB", 
                    "meanA", "meanB","log2FC", "pval", "extra")
  mat = as.data.frame(mat)
  mat$pval = as.numeric(mat$pval)
  mat$log2FC = as.numeric(mat$log2FC)
  mat$log10pval = -log10(mat$pval + 10 ^ -400)
  mat$diffexpressed = "NO"
  num = which((mat$log2FC > 1) & (mat$pval < 0.05))
  mat[num , "diffexpressed"] = "UP"
  
  num = which((mat$log2FC < -1) & (mat$pval < 0.05))
  mat[num , "diffexpressed"] = "DOWN"
  mat$diffexpressed = factor(mat$diffexpressed)
  
  png(paste0(clus1, " vs ", clus2, ".png"), width = 1000, height = 1000)
  plot(mat$log2FC, mat$log10pval, col = 
         c("red", "black", "blue")[mat$diffexpressed], pch = 20,
       xlab = "log2 foldchange", ylab = "-log10 pvalue", 
       main = paste0(clus1, " vs ", clus2))
  dev.off()
  
  write.csv(mat,file = paste0(clus1, " vs ", clus2, ".csv"))
}