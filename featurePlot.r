TransformCoords <- function(data) {  
  parentcoordinates = data[,1:2]
  childcoordinates = data[,3:4]  
  parentcentroid = c(mean(data[,1]), mean(data[,2]))
  childcentroid = c(mean(data[,3],na.rm = TRUE), mean(data[,4],na.rm = TRUE))  
  parentdist = sqrt((parentcoordinates[,1]-parentcentroid[1])^2 + (parentcoordinates[,2]-parentcentroid[2])^2)
  sparentdist = subset(parentdist, parentdist < quantile(parentdist,0.95))
  cutoff = mean(sparentdist)+3*sd(sparentdist)
  sparentcoordinates = parentcoordinates[which(parentdist < cutoff),]  
  bbparent = c(min(sparentcoordinates[,1], na.rm = TRUE), min(sparentcoordinates[,2], na.rm = TRUE),max(sparentcoordinates[,1], na.rm = TRUE) , max(sparentcoordinates[,2], na.rm = TRUE))  
  bbchild = c(min(childcoordinates[,1], na.rm = TRUE), min(childcoordinates[,2], na.rm = TRUE),max(childcoordinates[,1], na.rm = TRUE) , max(childcoordinates[,2], na.rm = TRUE))  
  wp = bbparent[3]-bbparent[1]
  hp = bbparent[4]-bbparent[2]
  wc = bbchild[3]-bbchild[1]
  hc = bbchild[4]-bbchild[2]  
  shrinkparent = 1.5
  sh = hc/(hp/shrinkparent)
  sw = wc/(wp/shrinkparent)  
  bbchildnew = childcoordinates
  bbchildnew[,1] = (childcoordinates[,1] - bbchild[1])/sw + bbparent[1]
  bbchildnew[,2] = (childcoordinates[,2] - bbchild[2])/sh + bbparent[2]  
  bbchildnew
}

#####################################################################################
library(Seurat)
library(ggplot2)
library(ggpubr)

featureSel = c("PTPRC", "CD14", "TMEM119", "CD163", "LYZ", "CD3E", "NKG7", "CD79A", "IFITM2", "S100A8", "PLP1", "PTPRZ1", "OLIG2", "MKI67", "CHI3L1", "VEGFA", "HILPDA", "PDGFRA", "EGFR", "DCN", "STMN2", "GFAP", "HAMP")

files = list.files(".",  pattern = "*.csv")
num = grep("UMAPLevel",  files)
files = files[num]

x <- sapply(files, function(x) strsplit(x, "_")[[1]][1], USE.NAMES=FALSE)
x = unique(x)

coords = matrix(NA, ncol = 2, nrow = length(x))

for(i in 2:length(x)){
  num = grep(x[i],  files)
  coords[i, 1] = min(num)
  coords[i, 2] = max(num)
}
coordsF = coords[-1,]

info = readRDS("OutputRCC1.rds")
#rownames(info) = info[,1]
#info = info[,-1]

pumap = read.csv("UMAPLevel1.csv")
rownames(pumap) = pumap[,1]
pumap = pumap[,-1]

data = readRDS("Level1.rds")

identical(rownames(data@meta.data),  rownames(info))
data@meta.data = info

identical(rownames(data@reductions$umap@cell.embeddings),  rownames(pumap))
data@reductions$numap = data@reductions$umap
data@reductions$numap@cell.embeddings = as.matrix(pumap)
colnames(data@reductions$numap@cell.embeddings) = c("UMAP_1",  "UMAP_2")

data@active.ident = factor(data@meta.data$Level_1)
names(data@active.ident) = rownames(data@meta.data)

meta = merge(as.matrix(data@meta.data), as.matrix(data@reductions$numap@cell.embeddings), by = 0, all = T)
meta = as.data.frame(meta)
lmarkers = list()

for(markers in 1:length(featureSel)){
  meta[, featureSel[markers]] = as.numeric(meta[, featureSel[markers]])
  print(class(meta[, featureSel[markers]]))
  p1 = ggplot(meta, aes_string(x = "UMAP_1", y = "UMAP_2", color = featureSel[markers])) + geom_point(size = 0.01) + 
    theme_classic() + theme(plot.title = element_text(size = 30)) + 
    scale_color_gradientn(colors = c("lightgrey", "red")) + ggtitle(featureSel[markers])
  
  #png(paste0("Level_1_", featureSel[markers], ".png"))
    #print(p1)
  #dev.off()
  lmarkers[[markers]] = p1  
}

png("Level1_markers.png", width = 3000, height = 2000)
  print(ggarrange(plotlist = lmarkers, ncol = 5, nrow = 5))
dev.off()

cols = "Level_1"
for(counts in 1:nrow(coordsF)){
  start = coordsF[counts, 1]
  end = coordsF[counts, 2]
  for(i in start:end){
    cumap = read.csv(files[i])
    rownames(cumap) = cumap[,1]
    cumap = cumap[,-1]
    
    c1umap = merge(pumap, cumap, by = 0, all.y = T)
    rownames(c1umap) = c1umap[,1]
    c1umap = c1umap[,-1]
    resUMAP = TransformCoords(c1umap)
    
    for(coords in 1:nrow(resUMAP)){
      pumap[as.character(rownames(resUMAP)[coords]), "parentUMAP_1"] = resUMAP[coords, 1]
      pumap[as.character(rownames(resUMAP)[coords]), "parentUMAP_2"] = resUMAP[coords, 2]
    }
  }
  data@reductions$numap@cell.embeddings = as.matrix(pumap)
  colnames(data@reductions$numap@cell.embeddings) = c("UMAP_1",  "UMAP_2")
  meta = data@meta.data
  cols = append(cols, paste0("Level_",  counts + 1))
  meta$newClus = apply(meta[ ,cols], 1, paste0, collapse = "_" )
  
  data@active.ident = factor(meta$newClus)
  names(data@active.ident) = rownames(meta)
  
  meta = merge(as.matrix(data@meta.data), as.matrix(data@reductions$numap@cell.embeddings), by = 0, all = T)
  meta = as.data.frame(meta)
  lmarkers = list()
  
  for(markers in 1:length(featureSel)){
    meta[, featureSel[markers]] = as.numeric(meta[, featureSel[markers]])
    print(class(meta[, featureSel[markers]]))
    p1 = ggplot(meta, aes_string(x = "UMAP_1", y = "UMAP_2", color = featureSel[markers])) + geom_point(size = 0.01) + 
      theme_classic() + theme(plot.title = element_text(size = 30)) + 
      scale_color_gradientn(colors = c("lightgrey", "red")) + ggtitle(featureSel[markers])
    png(paste0("Level_",  counts +  1, "_", featureSel[markers], ".png"))
      print(p1)
    dev.off()
    lmarkers[[markers]] = p1  
  }
  
  png(paste0("Level", counts + 1,"_markers.png"), width = 3000, height = 2000)
  print(ggarrange(plotlist = lmarkers, ncol = 5, nrow = 5))
  dev.off()
}



