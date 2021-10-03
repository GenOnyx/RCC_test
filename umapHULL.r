levels = unique(info$Level_1)

pdf("newUMAP_colored.pdf", width = 60, height = 60)

plot(info$UMAP_1, info$UMAP_2, pch = 20, cex = 0.1, 
     xlab = "UMAP_1", ylab = "UMAP_2", main = "Level_1")

levelcolors = c('blue', 'palevioletred4', 'orange', 
                'red', 'purple', 'hotpink', 'green', 'lightblue', 'brown',
                'darkgreen', 'darkslategray', 'deeppink3', 'plum2')

lines = c(1,2,3,4,5,6)

for(i in 1:length(levels)){
  i1 = subset(info, info$Level_1 == levels[i])
  #parentcentroid = c(mean(i1[,"UMAP_1"]), mean(i1[,"UMAP_2"]))
  
  #parentdist = sqrt((i1[,"UMAP_1"]-parentcentroid[1])^2 + (i1[,"UMAP_2"]-parentcentroid[2])^2)
  #sparentdist = subset(parentdist, parentdist < quantile(parentdist,0.95))
  #cutoff = mean(sparentdist)+3*sd(sparentdist)
  #i1 = i1[which(parentdist < cutoff),]
  
  X = as.matrix(i1[, c("UMAP_1", "UMAP_2")])
  a = chull(X)
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts,], col = levelcolors[levels[i]], lwd = 3, lty = 1)
}


#plot(info$UMAP_1, info$UMAP_2, pch = 20, cex = 0.1, xlab = "UMAP_1",
     #ylab = "UMAP_2", main = "Level_2")
levels = unique(info$l2)

for(i in 1:length(levels)){
  i1 = subset(info, info$l2 == levels[i])
  
  parentcentroid = c(mean(i1[,"UMAP_1"]), mean(i1[,"UMAP_2"]))
  
  parentdist = sqrt((i1[,"UMAP_1"]-parentcentroid[1])^2 + (i1[,"UMAP_2"]-parentcentroid[2])^2)
  sparentdist = subset(parentdist, parentdist < quantile(parentdist,0.95))
  cutoff = mean(sparentdist)+3*sd(sparentdist)
  i1 = i1[which(parentdist < cutoff),]
  
  lev1 = unique(i1$Level_2)
  if(lev1 == 0){ next}
  col = levelcolors[lev1]
  print(paste0(levels[i], " ", lev1, " ", col))
  
  X = as.matrix(i1[, c("UMAP_1", "UMAP_2")])
  a = chull(X)
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts,], col = col, lwd = 2, lty = 2)
}


plot(info$UMAP_1, info$UMAP_2, pch = 20, cex = 0.1, xlab = "UMAP_1",
     ylab = "UMAP_2", main = "Level_3")
levels = unique(info$l3)

for(i in 1:length(levels)){
  i1 = subset(info, info$l3 == levels[i])
  
  parentcentroid = c(mean(i1[,"UMAP_1"]), mean(i1[,"UMAP_2"]))
  
  parentdist = sqrt((i1[,"UMAP_1"]-parentcentroid[1])^2 + (i1[,"UMAP_2"]-parentcentroid[2])^2)
  sparentdist = subset(parentdist, parentdist < quantile(parentdist,0.95))
  cutoff = mean(sparentdist)+3*sd(sparentdist)
  i1 = i1[which(parentdist < cutoff),]
  
  lev1 = unique(i1$Level_3)
  if(lev1 == 0){ next}
  col = levelcolors[lev1]
  print(paste0(levels[i], " ", lev1, " ", col))
  
  X = as.matrix(i1[, c("UMAP_1", "UMAP_2")])
  a = chull(X)
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts,], col = col, lwd = 1.5, lty = 3)
}


plot(info$UMAP_1, info$UMAP_2, pch = 20, cex = 0.1, xlab = "UMAP_1",
     ylab = "UMAP_2", main = "Level_4")
levels = unique(info$l4)

for(i in 1:length(levels)){
  i1 = subset(info, info$l4 == levels[i])
  
  parentcentroid = c(mean(i1[,"UMAP_1"]), mean(i1[,"UMAP_2"]))
  
  parentdist = sqrt((i1[,"UMAP_1"]-parentcentroid[1])^2 + (i1[,"UMAP_2"]-parentcentroid[2])^2)
  sparentdist = subset(parentdist, parentdist < quantile(parentdist,0.95))
  cutoff = mean(sparentdist)+3*sd(sparentdist)
  #i1 = i1[which(parentdist < cutoff),]
  
  lev1 = unique(i1$Level_4)
  if(lev1 == 0){ next}
  col = levelcolors[lev1]
  print(paste0(levels[i], " ", lev1, " ", col))
  
  X = as.matrix(i1[, c("UMAP_1", "UMAP_2")])
  a = chull(X)
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts,], col = col, lwd = 1.5, lty = 4)
}


plot(info$UMAP_1, info$UMAP_2, pch = 20, cex = 0.1, xlab = "UMAP_1",
     ylab = "UMAP_2", main = "Level_5")
levels = unique(info$l5)

for(i in 1:length(levels)){
  i1 = subset(info, info$l5 == levels[i])
  
  parentcentroid = c(mean(i1[,"UMAP_1"]), mean(i1[,"UMAP_2"]))
  
  parentdist = sqrt((i1[,"UMAP_1"]-parentcentroid[1])^2 + (i1[,"UMAP_2"]-parentcentroid[2])^2)
  sparentdist = subset(parentdist, parentdist < quantile(parentdist,0.95))
  cutoff = mean(sparentdist)+3*sd(sparentdist)
  #i1 = i1[which(parentdist < cutoff),]
  
  lev1 = unique(i1$Level_5)
  if(lev1 == 0){ next}
  col = levelcolors[lev1]
  print(paste0(levels[i], " ", lev1, " ", col))
  
  X = as.matrix(i1[, c("UMAP_1", "UMAP_2")])
  a = chull(X)
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts,], col = col, lwd = 1.5, lty = 5)
}


plot(info$UMAP_1, info$UMAP_2, pch = 20, cex = 0.1, xlab = "UMAP_1",
     ylab = "UMAP_2", main = "Level_6")
levels = unique(info$l6)

for(i in 1:length(levels)){
  i1 = subset(info, info$l6 == levels[i])
  
  parentcentroid = c(mean(i1[,"UMAP_1"]), mean(i1[,"UMAP_2"]))
  
  parentdist = sqrt((i1[,"UMAP_1"]-parentcentroid[1])^2 + (i1[,"UMAP_2"]-parentcentroid[2])^2)
  sparentdist = subset(parentdist, parentdist < quantile(parentdist,0.95))
  cutoff = mean(sparentdist)+3*sd(sparentdist)
  #i1 = i1[which(parentdist < cutoff),]
  
  lev1 = unique(i1$Level_6)
  if(lev1 == 0){ next}
  col = levelcolors[lev1]
  print(paste0(levels[i], " ", lev1, " ", col))
  
  X = as.matrix(i1[, c("UMAP_1", "UMAP_2")])
  a = chull(X)
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts,], col = col, lwd = 1.5, lty = 6)
}
legend("topright", legend = c("K1", "K2", "K3", "K4", 
                             "K5", "K6", "K7", "K8", "K9", "K10", 
                             "K11", "K12", "K13"), col = levelcolors,
       cex=0.3, pch = 20, xpd = F)
legend("bottom", legend = paste0("Level_", 1:3), lty = c(1:3), xpd = F,
       cex = 0.3)

legend("topleft", legend = val, col = code, cex= 0.3, pch = 20)
plot(X, col = "white")
