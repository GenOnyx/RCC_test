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

########################################################################################################################################################################################################
arrange = function(data){
  obj = CreateSeuratObject(counts = data, project = "reference", min.cells = 3, min.features = 1)
  print(obj)

  dimMax = 10
  resMax = 0.1
  if(ncol(obj@assays$RNA@counts) < 150){
    dimMax = 5
  }
  if(ncol(obj@assays$RNA@counts) < 100){
    dimMax = 2
  }

  obj1 = obj
  med_obj = median(obj1@meta.data$nFeature_RNA)
  if(med_obj < 1000){
    summary1 = summary(obj1@meta.data$nFeature_RNA)[[5]]
    varGenes = round(summary1*0.2, 0)
  } else {
    varGenes = 1000
  }
  if(varGenes < 50){
    varGenes = 50
  }
  obj1 = NormalizeData(obj1, normalization.method = "LogNormalize", scale.factor = 100000)
  obj1 = FindVariableFeatures(obj1, selection.method = "vst", nfeatures = varGenes)
  all.genes = rownames(obj1)
  obj1 = ScaleData(obj1, features = all.genes)
  rIN = round(runif(8, min = 4, max = 12),0)
  obj1 = RunPCA(obj1, features = VariableFeatures(object = obj1), npcs = dimMax, seed.use = (runif(1, min = 1, max = 10^(rIN)))/10000)
  obj1 = FindNeighbors(obj1, dims = 1:dimMax)
  if(ncol(data) < 20000){
    algo = 4 
  } else {
    algo = 1
  }
  obj1 = FindClusters(obj1, resolution = 0.1, algorithm = algo, n.iter = 100)
  obj1 = RunUMAP(obj1, dims = 1:dimMax)
  meta = obj1@meta.data
  meta$seurat_clusters = as.numeric(meta$seurat_clusters)
  
  obj1@meta.data = meta
  obj1@active.ident = factor(obj1@meta.data$seurat_clusters)
  names(obj1@active.ident) = rownames(obj1@meta.data)
  checkOBJ <<- obj1 
  if(length(unique(obj1@meta.data$seurat_clusters)) > 1){
    markers = FindAllMarkers(obj1, features = obj1@assays$RNA@var.features, 
                             logfc.threshold = 1, min.pct = 0.50, only.pos = T)
    if(nrow(markers) == 0){
      obj1 = 0
      return(obj1)
    }
  }
  return(obj1)
}

library(Seurat)
library(ggplot2)
library(ggpubr)

RCC_clus = function(config_file){
  
  library(data.table)
  
  config_file = "config.csv"
  conf_file = read.csv(config_file, header = FALSE)
  if(grepl(".rds$", as.character(conf_file[1, 2]))){
    expr_data = readRDS(as.character(conf_file[1, 2]))
  } else {
    expr_data = as.data.frame(fread(as.character(conf_file[1, 2])))
    rownames(expr_data) = expr_data[, 1]
    expr_data = expr_data[, -1]
  }
  
  orig.data = expr_data
  print("matrix file read")
  
  sample_ids = colnames(expr_data)
  cluster_num = rep(1,ncol(expr_data))
  sample_ids = cbind(sample_ids,cluster_num)
  rownames(sample_ids) = sample_ids[,1]
  sample_ids = as.data.frame(sample_ids)
  
  m = NULL
  
  sampleInfo = as.data.frame(fread(as.character(conf_file[2,2])))
  rownames(sampleInfo) = sampleInfo[,1]
  sampleInfo = sampleInfo[,-1]
  init_cols = colnames(sampleInfo)
  init_colNum = ncol(sampleInfo)
  print("SampleInfo file read")
  
  common_samples = intersect(rownames(sample_ids),rownames(sampleInfo))
  sampleInfo = subset(sampleInfo,rownames(sampleInfo) %in% common_samples)
  sampleInfo_col_num = ncol(sampleInfo)
  sampleInfo <<- sampleInfo
  
  tcol_num = ncol(sampleInfo) + 1
  
  variableGenes <<- NULL
  min_samples = as.numeric(as.vector(conf_file[3,2]))
  finalOut <<- as.character(conf_file[4,2])
  
  colnames(sampleInfo) = init_cols
  
  if(dir.exists(finalOut)){
    print("Output directory found")
  } else {
    system(paste0('mkdir ', finalOut))
    print("Output directory created")
  }
  pdf(paste0(finalOut,"/UMAPres.pdf"), width = 20, height = 10)
  recursive_consensus = function(sample_ids,col_num){
   if((nrow(sample_ids) < min_samples)){
      selected_cluster = 1
    } else {
      expr_subset = NULL
      num = which(colnames(expr_data) %in% rownames(sample_ids))
      expr_subset = expr_data[, num]
      
      cluster_parallel <<- arrange(expr_subset)
      if(class(cluster_parallel) == "numeric"){
        selected_cluster = 0
      } else {
      variableGenes <<- append(variableGenes, cluster_parallel@assays$RNA@var.features)
      
      p1 = DimPlot(cluster_parallel, label = T, pt.size = 0.1)
      print(p1)
      cluster_file = cluster_parallel@meta.data
      cluster_file$seurat_clusters = as.numeric(cluster_file$seurat_clusters)
      cluster_file <<- cluster_file
      selected_cluster = max(cluster_file$seurat_clusters)
      }
    }
    print(paste0("Selected cluster is: ", selected_cluster))
    
    if(selected_cluster > 1){
      for(h in 1:nrow(cluster_file)){
        sample_name = as.character(rownames(cluster_file)[h])
        sampleInfo[sample_name,col_num] <<- cluster_file[sample_name,"seurat_clusters"]
      }
      sampleInfo <<- sampleInfo
      
      sampleLevel = abs(sampleInfo_col_num - col_num)
      if(sampleLevel == 1){
        saveRDS(cluster_parallel, file = paste0(finalOut, "/Level1.rds"))
        pumap <<- cluster_parallel@reductions$umap@cell.embeddings
        colnames(pumap) = c("parentUMAP_1", "parentUMAP_2")
        resUMAP = pumap
        cumap = pumap
      } else {
        cumap = cluster_parallel@reductions$umap@cell.embeddings
        c1umap = merge(pumap, cumap, by = 0, all.y = T)
        rownames(c1umap) = c1umap[,1]
        c1umap = c1umap[,-1]
        
        resUMAP = TransformCoords(c1umap)
      }
      
      for(coords in 1:nrow(resUMAP)){
        pumap[as.character(rownames(resUMAP)[coords]), "parentUMAP_1"] = resUMAP[coords, 1]
        pumap[as.character(rownames(resUMAP)[coords]), "parentUMAP_2"] = resUMAP[coords, 2]
      }
      pumap <<- pumap
      
      fileName = paste0("Level", sampleLevel)
      
      if(sampleLevel > 1){
        for(g in 1:(sampleLevel - 1)){
          sampleSel = as.character(rownames(cluster_file)[1])
          a = unique(sampleInfo[sampleSel, sampleInfo_col_num + g])
          fileName = paste(fileName,"_k",a, sep = "")
        }
      }
      library(ggplot2)
      library(ggpubr)
      
      m = col_num + 1
      
      write.csv(cumap, file = paste0(finalOut,"/UMAP", fileName, ".csv"))
      
      for(b in 1:selected_cluster){
        sample_ids = subset(cluster_file,cluster_file[,"seurat_clusters"] == b)
        sample_ids <<- sample_ids
        recursive_consensus(sample_ids,m)
      }
    }
  }
  recursive_consensus(sample_ids,tcol_num)
  
  number_of_levels = abs(sampleInfo_col_num - ncol(sampleInfo))
  print(number_of_levels)
  sampleInfo[,c((ncol(sampleInfo)-number_of_levels):ncol(sampleInfo))][is.na(sampleInfo[,c((ncol(sampleInfo)-number_of_levels):ncol(sampleInfo))])] = 0
  colnames(sampleInfo)[(sampleInfo_col_num + 1): ncol(sampleInfo)] = paste0("Level_",1:number_of_levels)
  colNum = c(1:number_of_levels)
  cols = paste("Level_", colNum, sep = "")
  
  if(length(cols) == 1){
    Clusters = sampleInfo[,ncol(sampleInfo)]
    sampleInfo = cbind(sampleInfo,Clusters)
    colnames(sampleInfo)[ncol(sampleInfo)] = "Clusters"
  } else {
    colnames(sampleInfo)[c(((ncol(sampleInfo) - number_of_levels)+1):ncol(sampleInfo))] = cols
    sampleInfo[,(ncol(sampleInfo) + 1)] = do.call(paste0, sampleInfo[c(cols)])
    colnames(sampleInfo)[ncol(sampleInfo)] = "Concatenated_clusters"
    
    len = length(unique(sampleInfo[,ncol(sampleInfo)]))
    cluster_found = unique(sampleInfo[,ncol(sampleInfo)])
    
    for(k in 1:len){
      for(j in 1:nrow(sampleInfo)){
        if(sampleInfo[j,ncol(sampleInfo)] == cluster_found[k]){
          sampleInfo[j,ncol(sampleInfo)] = k
        }
      }
    }
    
    colnames(sampleInfo)[ncol(sampleInfo)] = "Clusters"
    sampleInfo <<- sampleInfo
  }
  samples = intersect(colnames(orig.data), rownames(sampleInfo))
  num = which(colnames(orig.data) %in% samples)
  orig.data = orig.data[, num]
  
  sampleInfo = subset(sampleInfo, rownames(sampleInfo) %in% samples)
  orig.data = orig.data[, order(colnames(orig.data))]
  sampleInfo = sampleInfo[order(rownames(sampleInfo)), ]	
  
  obj = CreateSeuratObject(counts = orig.data, project = "RCCprocessed", min.cells = 3, min.features = 1, meta.data = sampleInfo)
  saveRDS(obj, file = paste0(finalOut, "/OutputRCC.rds"))
  write.csv(sampleInfo, file = paste0(finalOut, "/OutputRCC.csv"))
  write.csv(pumap,file = paste0(finalOut, "/UMAPcoords.csv"))
  dev.off()

  write.csv(unique(variableGenes), file = paste0(finalOut, "/genesUsed.csv"))
  source("/media/user/disk3/R/clusterAnnotation.r")
  clusterAnnotation("OutputRCC1.csv",  "ID", "Clusters", finalOut)
}