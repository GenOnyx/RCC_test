library(Matrix)
library(rhdf5)
library(tidyverse)
library(glue)
library(Seurat)

retCounts = function(h5file){
  h5 = h5read(h5file, "matrix")
  counts = sparseMatrix(
    dims = h5$shape,
    i = as.numeric(h5$indices),
    p = as.numeric(h5$indptr),
    x = as.numeric(h5$data),
    index1 = FALSE
  )
  colnames(counts) = h5$barcodes
  rownames(counts) = h5$features$name
  return(counts)
}

files = list.files(".", pattern = "*.h5")

for(i in 1:length(files)){
  print(i)
  counts = retCounts(files[i])
  obj =  CreateSeuratObject(counts = counts, colData = info, project = "RNAseq", min.cells = 3, min.features = 50)
  genes = c("IFITM2", "S100A12", "S100A8")
  counts = as.matrix(obj@assays$RNA@counts)
  m = subset(counts, rownames(counts) %in% genes)
  dim(m)
  
  IFITM2_S100A12_S100A8 = data.frame(apply(m, 2, mean))
  colnames(IFITM2_S100A12_S100A8) = "IFITM2_S100A12_S100A8"
  mat = obj@meta.data
  mat = cbind.data.frame(mat, IFITM2_S100A12_S100A8)
  
  Select = rep(0, nrow(mat))
  num = which((mat$nCount_RNA > 1000) | (mat$IFITM2_S100A12_S100A8 > 2))
  Select[num] = 1
  mat = cbind.data.frame(mat, Select)
  obj@meta.data = mat
  obj1 = subset(x = obj, subset = Select == 1)
  
  obj1
  #obj1 = NormalizeData(obj1, normalization.method = "LogNormalize", scale.factor = 10000)
  #obj1 = FindVariableFeatures(obj1, selection.method = "vst", nfeatures = 2000)
  #all.genes = rownames(obj1)
  #obj1 = ScaleData(obj1, features = all.genes)
  #obj1 = RunPCA(obj1, features = VariableFeatures(object = obj1))
  #obj1 = FindNeighbors(obj1, dims = 1:10)
  #obj1 = FindClusters(obj1, resolution = 0.5)
  #obj1 = RunUMAP(obj1, dims = 1:10)
  #DimPlot(obj1, reduction = "umap", label = T)
  
  saveRDS(obj1, file = gsub(".h5", ".rds",files[i]))
  #genes = c("CD45", "CD14", "TMEM119", "CD163", "LYZ", "CD3E", "NKG7", "CD79A", "IFITM2", "S100A8", "PLP1", "PTPRZ1", "OLIG2", "MKI67", "CHI3L1", "VEGFA", "HILPDA", "PDGFRA", "EGFR", "DCN", "STMN2", "GFAP", "HAMP", "HBB", "PTPRC", "VSNL1", "DCX")
  
  #pdf(gsub(".h5", ".pdf",files[i]), height = 30, width = 20)
  #	FeaturePlot(obj1, features = genes, cols = c("grey", "red"))
  #dev.off()
}

files = list.files(".", pattern = "*.rds")

mat_BSC = NULL
info_BSC = NULL
mat_core = NULL
info_core = NULL
mat_peri = NULL
info_peri = NULL

for(i in 1:length(files)){
  if(grepl("BSC", files[i])){
    obj1 = readRDS(files[i])
    counts = as.matrix(obj1@assays$RNA@counts)
    colnames(counts) = paste0(gsub(".rds", "", files[i]), "_", colnames(counts))
    info = obj1@meta.data
    rownames(info) = paste0(gsub(".rds", "", files[i]), "_", rownames(info))
    if(is.null(mat_BSC)){
      mat_BSC = counts
    } else {
      mat_BSC = merge(mat_BSC, counts, all = T, by = 0)
      rownames(mat_BSC) = mat_BSC[,1]
      mat_BSC = mat_BSC[,-1]
    }
    info_BSC = rbind(info_BSC, info)
  }
  if(grepl("TSC", files[i])){
    if(grepl("Core", files[i])){
      obj1 = readRDS(files[i])
      counts = as.matrix(obj1@assays$RNA@counts)
      colnames(counts) = paste0(gsub(".rds", "", files[i]), "_", colnames(counts))
      info = obj1@meta.data
      rownames(info) = paste0(gsub(".rds", "", files[i]), "_", rownames(info))
      if(is.null(mat_core)){
        mat_core = counts
      } else {
        mat_core = merge(mat_core, counts, all = T, by = 0)
        rownames(mat_core) = mat_core[,1]
        mat_core = mat_core[,-1]
      }
      info_core = rbind(info_core, info)
    }
    if(grepl("Peri", files[i])){
      obj1 = readRDS(files[i])
      counts = as.matrix(obj1@assays$RNA@counts)
      colnames(counts) = paste0(gsub(".rds", "", files[i]), "_", colnames(counts))
      info = obj1@meta.data
      rownames(info) = paste0(gsub(".rds", "", files[i]), "_", rownames(info))
      if(is.null(mat_peri)){
        mat_peri = counts
      } else {
        mat_peri = merge(mat_peri, counts, all = T, by = 0)
        rownames(mat_peri) = mat_peri[,1]
        mat_peri = mat_peri[,-1]
      }
      info_peri = rbind(info_peri, info)
    }
    
  }
}
