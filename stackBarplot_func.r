stackBarplot = function(infoFile, attr1, attr2, colorFile){
  library(RColorBrewer)
  info = read.csv(infoFile)
  rownames(info) = info[,1]
  info = info[,-1]
  
  a = NULL
  cellTypes = unique(info[,  attr1])
  
  for(i in 1:length(cellTypes)){
    i1 = subset(info, info[, attr1] == cellTypes[i]) 
    b = as.data.frame(table(i1[, attr2]))
    b = cbind.data.frame(rep(cellTypes[i], nrow(b)), b)
    b[,4] = (b[,3]/sum(b[,3])) * 100
    colnames(b) = c("RCC_clus", "ID", "Value", "Percent")
    a = rbind(a, b)
  }
 # a[,2] = factor(a[,2])
  if(colorFile == "none"){
    colors = rainbow(length(unique(a[,2])))
    print("no custom colors used")
  } else {
    colors = read.csv(colorFile)
    colors = subset(colors, colors[,1] %in% a[,2])
    colors = colors[,2]
  }
  #b = merge(a, colors, by = "ID", all.x = T)
  
  if(attr2 == "sample_location"){
    png(paste0(attr1, "_", attr2, "barplot.png"),width = 2000,  height = 1000)
    print(ggplot(a, aes(fill = ID, y = Percent, x = RCC_clus)) + geom_bar(position = "stack", stat = "identity") + 
            theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1)) + 
            scale_fill_manual(values = c("Core_SNU27_Core" = "lightblue",
                                "Core_SNU33_Core" = "pink",
                                "Core_SNU33b2_Core" = "hotpink",
                                "Core_SNU40_Core" = "springgreen",
                                "Core_SNU43_Core" = "yellow",
                                "Accutase_UK" = "purple",
                                "GBM_v7_UK" = "mediumorchid4",
                                "Peri_SNU27_Peri" = "navy",
                                "Peri_SNU33_Peri" = "hotpink3",
                                "Peri_SNU33b2_Peri" = "hotpink4",
                                "Peri_SNU40_Peri" = "olivedrab",
                                "Peri_SNU43_Peri" = "goldenrod1",
                                "SNU38_Control" = "black")))
    dev.off()
  } else {
    png(paste0(attr1, "_", attr2, "barplot.png"),width = 2000,  height = 1000)
    print(ggplot(a, aes(fill = ID, y = Percent, x = RCC_clus)) + geom_bar(position = "stack", stat = "identity") + 
            theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1)) + 
            scale_fill_manual(values = colors))
    dev.off()
  }
}
