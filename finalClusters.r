library(gsubfn)
library(matrixStats)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(RecursiveConsensusClustering)

matrixAll = "../count.rds"
trackingPlot(matrixAll, "OutputRCC.csv", "genesUsed.csv", "attr", ".")
fileName = matrixAll


files = list.files(".", pattern = "Level_")
files = files[!grepl("*k_[0-9].csv", files)]
files = files[!grepl("*k_10.csv", files)]
files = sort(str_replace(files, "10", "T"))
final_mat = NULL
topGenes = 10

for(i in 1:length(files)){
	mat = NULL
	
	files[i] = str_replace(files[i], "T", "10")
	
	data = read.csv(files[i])
	data = data[,-1]
	for(j in 1:ncol(data)){
		mat = append(mat, as.character(data[c(1:topGenes), j]))
	}
	files1 = str_replace(files[i], "Level_", "")
	files2 = str_replace(files1, ".csv", "")
	files3= str_replace(files2, "k1", "")
	lev = strsplit(files3, "_")[[1]][1]
	levAll = rep(lev, ncol(data) * topGenes)
	if(lev == "1"){
		colnames(data) = str_replace(colnames(data), "X", "")
		clus = rep(colnames(data), each = topGenes)
	} else {
		prevClus = strsplit(files3, "_")[[1]][2]
		colnames(data) = str_replace(colnames(data), "X", "")
		clus = rep(paste0(prevClus, colnames(data)), each = topGenes)
	}
	final = cbind(mat, levAll)
	final = cbind(final, clus)
	final_mat = rbind(final_mat, final)

}
final_mat[final_mat == "0"] <- NA
final_mat = final_mat[-which(is.na(final_mat[,1])),]
#final_mat = final_mat[-which(final_mat[,1] == "0"),]

colnames(final_mat) = c("genes", "Level", "clusters")
write.csv(final_mat, file = "finalClusters.csv", row.names = F)

if(grepl(".rds$", fileName)){
  data = readRDS(fileName)
} else{
  data = as.data.frame(fread(fileName))
  rownames(data) = data[,1]
  data = data[,-1]
}
sample = read.csv("sampleOrder.csv")
rownames(sample) = sample[,1]
sample = sample[,-1]

Levels = which(grepl("Level_*", colnames(sample)))
mLev = Levels[length(Levels)] + 1
sample = sample[, c(Levels[1]:mLev)]

for(i in 2:ncol(sample)){
	sample[,i] = str_replace(sample[,i], "10", "T")
}

gene = read.csv("finalClusters.csv")
a = paste0(gene[,1] , "_", gene[,3])
rownames(gene) = a


m = as.numeric(gene[,3])
mn = unique(gene[,3])

for(fk in 1:length(mn)){
	for(jk in 1:length(m)){
		if(m[jk] == mn[fk]){
			m[jk] = fk
		}
	}
}

gene = cbind.data.frame(gene, m)
gene1 = subset(gene, gene[,2] == 1)
gene2 = subset(gene, gene[,2] != 1)
gene2[,3] =  str_replace(gene2[,3], "10", "T")

gene = rbind(gene1, gene2)



mat = NULL
subData= subset(data, rownames(data) %in% as.character(gene[,1]))

for(i in 1:nrow(sample)){
	id = as.character(rownames(sample)[i])
	mat = cbind(mat, subData[,id])
}
colnames(mat) = rownames(sample)
rownames(mat) = rownames(subData)

fmat = NULL
for(i in 1:nrow(gene)){
	geneID = as.character(gene[i,1])
	fmat = rbind(fmat, mat[geneID, ])
}	

rownames(fmat) = rownames(gene)
colnames(fmat) = colnames(mat)

write.csv(fmat, file = "orderedData.csv")

col = apply(sample, 2, as.character)
rw = gene[,-1]
rw$Level = as.character(rw$Level)
rw$clusters = as.character(rw$clusters)
ha =  HeatmapAnnotation(df = col)
ha1 = rowAnnotation(df = rw)

zscore_mat = (fmat - rowMeans(fmat))/(rowSds(as.matrix(fmat)))[row(fmat)]
#ht = Heatmap(zscore_mat, col = colorRamp2(c(-2,0,2), c("orangered", "white", "purple")), top_annotation=ha, cluster_rows=F, cluster_columns=F)
#
#pdf("HeatmapMarker.pdf", height = 10, width = 15)
#	ha1 + ht
#dev.off()
#
#
#lev = unique(gene$Level)
#
mat = zscore_mat

#averaging for smooth plot

#for(i in 1:length(lev)){
#	print(paste0("Level_", i))
#	subGene = subset(gene, gene$Level == as.character(lev[i]))
#	subData = subset(zscore_mat, rownames(zscore_mat) %in% rownames(subGene))
#	subData = t(subData)
#	a = apply(as.data.frame(sample[,c(1:i)]), 1, paste, collapse = "_")
#	subData = cbind.data.frame(a, subData)
#	m = subData
#	mm = aggregate(m[,2:ncol(m)], list(m[,1]), median)
#	rownames(mm) = mm[,1]
#	mm = mm[,-1]
#	
#	df = as.data.frame(rownames(sample))
#	df = cbind.data.frame(df, a)
#	rownames(df) = df[,1]
#	
#	for(j in 1:nrow(mat)){
#		#print(paste0("row: ", j))
#		geneID = as.character(rownames(mat)[j])
#		if(geneID %in% colnames(mm)){
#			for(k in 1:ncol(mat)){
#				#print(paste0("col: ", k))
##what are we doing here??This loop takes forever
#				sampleID = as.character(colnames(mat)[k])
#				clus = df[sampleID, 2]
#				mat[j,k] = mm[as.character(clus), geneID]
#			}
#		}
#	}
#}

write.csv(mat, file = "zscoreAverage.csv")

pdf("HeatmapMarker_highres.pdf", height = 40, width = 60)
	col = apply(sample, 2, as.character)
	colorsR = c("black", "cyan", "green", "hotpink", "maroon", "yellow")
	Heatmap(mat, col = colorRamp2(c(-2,0,2), c("orangered", "white", "purple")), cluster_rows=F, cluster_columns=F, row_names_gp = gpar(fontsize = 0), column_names_gp = gpar(fontsize = 0), name = "ht")
	tabM = as.numeric(col[,3])
	tabM = as.data.frame(table(tabM))
	tabM = tabM[order(tabM[,1]),]
	tabM[,3] = (tabM[,2]/sum(tabM[,2]))
	tabM[,4] = cumsum(tabM[,3])
	
	tabN = as.numeric(rw[,1])
	tabN = as.data.frame(table(tabN))
	tabN = tabN[order(tabN[,1]),]
	tabN[,3] = (tabN[,2]/sum(tabN[,2]))
	tabN[,4] = cumsum(tabN[,3])
	tabN[,5] = 1 - tabN[,4]
	prev = 1
	tabK = as.numeric(rw[,3])
	tabK = as.data.frame(table(tabK))
	tabK = tabK[order(tabK[,1]),]
	tabK[,3] = (tabK[,2]/sum(tabK[,2]))
	tabK[,4] = cumsum(tabK[,3])
	tabK[,5] = 1 - tabK[,4]
	tabK = cbind(unique(rw[,2]), tabK)
	rownames(tabK) = tabK[,1]
	rownames(tabM) = tabM[,1]
	rownames(tabN) = tabN[,1]
	
	xPoint = NULL
	lev = unique(gene$Level)
	newCol = NULL
	perSample = 1/ncol(mat)
	col = cbind(col, 1:nrow(col))
	prevRow = 1
	
	for(row in 1:length(lev)){
		#mCol = colorsR[row]
		newCol = paste0(newCol, col[,row])
		col = cbind(col, newCol)
		newSub = subset(col, col[,row] != 0)
		nextLine = as.numeric(newSub[1, 1]) - 1
		prevLine = as.numeric(newSub[1, ncol(newSub)]) - 1
		prevClus = newSub[1,3]
		for(vert in 1:nrow(newSub)){
			if(newSub[vert, ncol(newSub)] != prevLine){ 
				print(paste0("Level: ", row, " prevRow: ", (as.numeric(newSub[vert,ncol(newSub) - 1]) - 1) * perSample, " Cluster: ", newSub[vert, ncol(newSub)]))
				print(decorate_heatmap_body("ht", { grid.lines((as.numeric(newSub[vert,ncol(newSub) - 1]) - 1) * perSample, c(as.numeric(tabN[row,5]), prevRow), gp = gpar(lty = 2, lwd = 1, col = "black"))}))
				prevClus = newSub[vert,3]
				if(is.null(xPoint)){
					xPoint = append(xPoint, (as.numeric(newSub[vert,ncol(newSub) - 1]) - 1) * perSample) 
				} else if((as.numeric(newSub[vert,ncol(newSub) - 1]) - 1) * perSample != xPoint[length(xPoint)]){
					
					xPoint = append(xPoint, (as.numeric(newSub[vert,ncol(newSub) - 1]) - 1) * perSample)	
					print(paste0("point added_1: ", xPoint[length(xPoint)]))
				}
				
			}
			if(vert == 1){
				prevLine = newSub[vert, ncol(newSub)]
				next
			}
			if(abs(as.numeric(newSub[(vert), 1]) - as.numeric(nextLine)) > 1){
				print(paste0("condition met ", (as.numeric(newSub[vert-1,ncol(newSub) - 1]) ) * perSample))
				print(decorate_heatmap_body("ht", { grid.lines((as.numeric(newSub[vert-1,ncol(newSub) - 1]) ) * perSample, c(as.numeric(tabN[row,5]), prevRow), gp = gpar(lty = 2, lwd = 1, col = "black"))}))	
				prevLine = newSub[vert, ncol(newSub)]
				prevClus = newSub[vert,3]
				
				print(paste0("last element: ", xPoint[length(xPoint)]))
				print(paste0("current xPoint: ", (as.numeric(newSub[vert-1,ncol(newSub) - 1]) - 1) * perSample))
				if(is.null(xPoint)){
					xPoint = append(xPoint, (as.numeric(newSub[vert-1,ncol(newSub) - 1]) - 1) * perSample)
				} else if((as.numeric(newSub[vert-1,ncol(newSub) - 1]) - 1) * perSample != xPoint[length(xPoint)]){
					
					xPoint = append(xPoint, (as.numeric(newSub[vert-1,ncol(newSub) - 1]) - 1) * perSample)	
					tmp = xPoint[length(xPoint)]
					xPoint[length(xPoint)] = xPoint[length(xPoint)-1]
					xPoint[length(xPoint)-1] = tmp
					print(paste0("point added_2: ", xPoint[length(xPoint)]))
				}
			}
			prevLine = newSub[vert, ncol(newSub)]
			nextLine = newSub[vert, 1]
		}
		print(paste0("Level: ", row, " prevRow: ", (as.numeric(newSub[vert,ncol(newSub) - 1])) * perSample, " Cluster: ", newSub[vert, ncol(newSub)]))
		print(decorate_heatmap_body("ht", { grid.lines((as.numeric(newSub[vert,ncol(newSub) - 1])) * perSample, c(as.numeric(tabN[row,5]), prevRow), gp = gpar(lty = 2, lwd = 1, col = "black"))}))
		
		xPoint = append(xPoint, (as.numeric(newSub[vert,ncol(newSub) - 1])) * perSample)
		print(paste0("point added_3: ", xPoint[length(xPoint)]))
		prevRow = as.numeric(tabN[row,5])
		col = col[, -c(ncol(col))]
		print(dim(col))
	}
	
	
	col = apply(sample, 2, as.character)
	y1 = 1
	val = 1
	lev = unique(gene$Level)
	newCol = NULL
	perSample = 1/ncol(mat)
	prevRow = 1
	
	index = c(1:nrow(col))
	col = cbind(col, index)
	Levels = which(grepl("Level_*", colnames(col)))
	a = col[, Levels[1]]
	col = cbind(col, a)
	
	col = apply(col, 2, function(x)gsub('\\s+', '', x))
	
	for(i in 2:length(Levels)){
		a = apply(col[, c(Levels[1]:Levels[i])], 1, paste0, collapse = "")
		col = cbind(col, a)
	}
	colnames(col)[c((ncol(col) - (length(Levels) -1)): ncol(col))] = paste0("newLevel_", Levels)
	Levels = which(grepl("newLevel_*", colnames(col)))
	kmn = 1
	
	for(knb in 1:length(Levels)){
		col[, Levels[knb]] = str_replace(col[, Levels[knb]], "10", "T")
	}
	LevelNumR = 1
	
	for(row in Levels[1]:Levels[length(Levels)]){
		mCol = colorsR[kmn]
		print(paste0("Value: ", val))
		clustersN = unique(col[,row])
		geneSub = subset(rw, rw[,1] == LevelNumR)
		for(clus in 1:length(clustersN)){
			newSub = subset(col, col[,row] == clustersN[clus])
			x1 = (as.numeric(newSub[1, "index"]) - 1) * perSample
			x2 = as.numeric(newSub[nrow(newSub),"index"]) * perSample
			clusNum = clustersN[clus]
			
			
			if((clusNum %in% rownames(tabK)) & (clusNum %in% geneSub[,2])){
				ind = which(rownames(tabK) == clusNum)
				if(ind == 1){
					y1 = 1
				} else {
					y1 = tabK[(ind - 1), 6]
				}
				y2 = tabK[clusNum, 6]
				h = abs(y1 - y2)
				w = abs(x1 - x2)
				gpx = (x2 - w/2)
				gpy = (y1 - h/2)
				print(paste0("clusterNum: ", clusNum, " x1: ", x1, " x2: ", x2, " y1: ", y1, " y2: ", y2))
				print(decorate_heatmap_body("ht", { grid.rect(x = gpx, y = gpy, height = h, width = w, gp = gpar(lwd = 4, col = mCol))}))
			}
		}
		kmn = kmn + 1
		LevelNumR = LevelNumR + 1
	}
dev.off()

source("/media/DSRG4new/RCC/attributeMarkers.r")
attributeMarkers(matrixAll, "OutputRCC.csv", "attr", "Clusters", 10)

