info = read.csv("OutputRCC.csv")
rownames(info) = info[,1]
info = info[,-1]

info$attr = paste0(info$sampleID, "_", info$location)

a = NULL

cellTypes = unique(info$cellType_l1)
for(i in 1:length(cellTypes)){
 i1 = subset(info, info$cellType_l1 == cellTypes[i]) 
 b = as.data.frame(table(i1$attr))
 b = cbind.data.frame(rep(cellTypes[i], nrow(b)), b)
 b[,4] = (b[,3]/sum(b[,3])) * 100
 colnames(b) = c("RCC_clus", "ID", "Value", "Percent")
 a = rbind(a, b)
}

ggplot(a, aes(fill = ID, y = Percent, x = RCC_clus)) + geom_bar(position = "stack", stat = "identity") + 
  theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  scale_fill_manual(values = c("red", "seagreen", "cyan", "yellow", "pink", "lightblue", "brown", "orange",
                               "black", "purple", "grey", "springgreen1", "violetred2"))


info$newClus = paste0(info$cellType_l1, "_", info$Level_2, "_", info$Level_3,  "_", info$Level_4, "_", 
                      info$Level_5, "_", info$Level_6, "_", info$Level_7)

a = NULL

cellTypes = unique(info$newClus)
for(i in 1:length(cellTypes)){
  i1 = subset(info, info$newClus == cellTypes[i]) 
  b = as.data.frame(table(i1$attr))
  b = cbind.data.frame(rep(cellTypes[i], nrow(b)), b)
  b[,4] = (b[,3]/sum(b[,3])) * 100
  colnames(b) = c("RCC_clus", "ID", "Value", "Percent")
  a = rbind(a, b)
}

ggplot(a, aes(fill = ID, y = Percent, x = RCC_clus)) + geom_bar(position = "stack", stat = "identity") + 
  theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  scale_fill_manual(values = c("red", "seagreen", "cyan", "yellow", "pink", "lightblue", "brown", "orange",
                               "black", "purple", "grey", "springgreen1", "violetred2"))
