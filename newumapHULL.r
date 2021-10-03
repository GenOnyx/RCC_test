data = readRDS("finalSeuratObj.rds")

info = read.csv("OutputRCC1.csv")
rownames(info) = info[,1]
info = info[,-1]

markers = read.csv("finalClusters.csv")

levels = unique(info$Level_1)

pdf("newUMAP_colored1.pdf", width = 60, height = 60)

plot(info$UMAP_1, info$UMAP_2, pch = 20, cex = 0.1, 
     xlab = "UMAP_1", ylab = "UMAP_2", main = "Level_1")

levelcolors = c('blue', 'palevioletred4', 'orange', 
                'red', 'purple', 'hotpink', 'green', 'lightblue', 'brown',
                'darkgreen', 'darkslategray', 'deeppink3', 'plum2')

lines = c(1,2,3,4,5,6)

for(i in 1:length(levels)){
  i1 = subset(info, info$Level_1 == levels[i])

  X = as.matrix(i1[, c("UMAP_1", "UMAP_2")])
  a = chull(X)
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts,], col = levelcolors[levels[i]], lwd = 3, lty = 1)
}

cinfo = info[, c("Level_1", "UMAP_1", "UMAP_2")]
minfo = aggregate(cinfo[, c(2,3)], list(cinfo[,1]), max)
ainfo = aggregate(cinfo[, c(2,3)], list(cinfo[,1]), mean)

finfo = cbind.data.frame(ainfo[, 1:2], minfo[,3])
colnames(finfo) = c("Clusters", "UMAP_1", "UMAP_2")
finfo[, 2:3] = finfo[, 2:3] + 0.2
with(finfo, text(UMAP_2 ~ UMAP_1, labels = Clusters, cex = 5, 
                 col = "grey50"))


levels = unique(info$l2)
finfo = NULL

for(i in 1:length(levels)){
  i1 = subset(info, info$l2 == levels[i])
  
  parentcentroid = c(mean(i1[,"UMAP_1"]), mean(i1[,"UMAP_2"]))
  
  parentdist = sqrt((i1[,"UMAP_1"]-parentcentroid[1])^2 + (i1[,"UMAP_2"]-parentcentroid[2])^2)
  sparentdist = subset(parentdist, parentdist < quantile(parentdist,0.95))
  cutoff = mean(sparentdist)+ 3 * sd(sparentdist)
  i1 = i1[which(parentdist < cutoff),]
  
  u2max = max(i1$UMAP_2)
  u1mean = mean(i1$UMAP_1)
  clus = unique(i1$Level_2)
  rowq = c(clus, u1mean, u2max)
  
  finfo = rbind(finfo, rowq)
  
  lev1 = unique(i1$Level_2)
  if(lev1 == 0){ next}
  col = levelcolors[lev1]
  print(paste0(levels[i], " ", lev1, " ", col))
  
  X = as.matrix(i1[, c("UMAP_1", "UMAP_2")])
  a = chull(X)
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts,], col = col, lwd = 1.5, lty = 2)
}

colnames(finfo) = c("labels", "UMAP_1", "UMAP_2")
finfo[, 2:3] = finfo[, 2:3] + 0.1
finfo = subset(finfo, finfo[,1] != 0)
finfo = as.data.frame(finfo)
with(finfo, text(UMAP_2 ~ UMAP_1, labels = labels, cex = 2, 
                 col = "grey50"))


levels = unique(info$l3)
finfo = NULL

for(i in 1:length(levels)){
  i1 = subset(info, info$l3 == levels[i])
  
  parentcentroid = c(mean(i1[,"UMAP_1"]), mean(i1[,"UMAP_2"]))
  
  parentdist = sqrt((i1[,"UMAP_1"]-parentcentroid[1])^2 + (i1[,"UMAP_2"]-parentcentroid[2])^2)
  sparentdist = subset(parentdist, parentdist < quantile(parentdist,0.95))
  cutoff = mean(sparentdist)+ 3 * sd(sparentdist)
  i1 = i1[which(parentdist < cutoff),]
  
  u2max = max(i1$UMAP_2)
  u1mean = mean(i1$UMAP_1)
  clus = unique(i1$Level_3)
  rowq = c(clus, u1mean, u2max)
  
  finfo = rbind(finfo, rowq)
  
  lev1 = unique(i1$Level_3)
  if(lev1 == 0){ next}
  col = levelcolors[lev1]
  print(paste0(levels[i], " ", lev1, " ", col))
  
  X = as.matrix(i1[, c("UMAP_1", "UMAP_2")])
  a = chull(X)
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts,], col = col, lwd = 0.75, lty = 5)
}

colnames(finfo) = c("labels", "UMAP_1", "UMAP_2")
finfo[, 2:3] = finfo[, 2:3] + 0.01
finfo = subset(finfo, finfo[,1] != 0)
finfo = as.data.frame(finfo)
with(finfo, text(UMAP_2 ~ UMAP_1, labels = labels, cex = 1, 
                 col = "grey50"))



levels = unique(info$colorLevel)

m1 = min(info$UMAP_1)
m2 = max(info$UMAP_1)

m3 = min(info$UMAP_2)
m4 = max(info$UMAP_2)

for(i in 1:length(levels)){
  i1 = subset(info, info$colorLevel == levels[i])
  col = unique(i1$colors)
  if(col == "black"){next}
  print(col)
  par(new = TRUE)
  plot(i1$UMAP_1, i1$UMAP_2, pch = 20, cex = 0.1, 
       xlab = "UMAP_1", ylab = "UMAP_2", main = "Level_1", xlim = c(m1,m2),
       ylim = c(m3,m4), col = col)
 
}
dev.off()

