rgb_plot = function(meta, m, nLev, fclus){
  library(magick)
  
  nLev = gsub("Level_", "", nLev)
  nLev = as.numeric(nLev)

  frink <- image_read("legend.png")
  
  max1 = max(meta$UMAP_1)
  min1 = min(meta$UMAP_1)
  
  max2 = max(meta$UMAP_2)
  min2 = min(meta$UMAP_2)
  
  diff1 = max1 - min1
  diff2 = max2 - min2
  
  
  if(nLev == 1){
    cex = 0.1
  } 
  if(nLev == 2){
    cex = 0.5
  }
  if(nLev == 3){
    cex = 1
  }
  if(nLev == 4){
    cex = 2
  }
  if(nLev == 5){
    cex = 3
  }
  if(nLev == 6){
    cex = 5
  }
  if(length(m) == 1){
    p1 = ggplot(meta, aes_string(x = "UMAP_1", y = "UMAP_2", color = m[1])) + 
      geom_point(size = 0.01)  + 
      theme(plot.title = element_text(size = 30)) + 
      scale_color_gradientn(colors = c("black", "red"))
    h = layer_data(p1)
    par(bg="gray40")
    plot(h$x, h$y, pch = 20,cex = cex, cex.lab = 5,
         ylab = "", col = h$colour, xlab = fclus, cex.main = 5,
         main = paste0("Level ", nLev, " Markers ", m[1]), 
         xlim = c(min1, max1 + 3))
    rasterImage(frink, max1, min2, max1 + (diff1/5), min2 + (diff2/5))
  } else if(length(m) == 2){
    p1 = ggplot(meta, aes_string(x = "UMAP_1", y = "UMAP_2", color = m[1])) + 
      geom_point(size = 0.01)  + 
      theme(plot.title = element_text(size = 30)) + 
      scale_color_gradientn(colors = c("black", "white"))
    
    p2 = ggplot(meta, aes_string(x = "UMAP_1", y = "UMAP_2", color = m[2])) + 
      geom_point(size = 0.01) + 
      theme(plot.title = element_text(size = 30)) + 
      scale_color_gradientn(colors = c("black", "white"))
    
    cols = data.frame(layer_data(p1)$colour,
                      layer_data(p2)$colour,
                      rep("#000000", nrow(meta)))
    
    cols$fcolor = NA
    for(i in 1:nrow(cols)){
      r = unlist(strsplit(cols[i,1], split = ""))[2:3]
      r = paste0(r[1],r[2])
      g = unlist(strsplit(cols[i,2], split = ""))[4:5]
      g = paste0(g[1],g[2])
      b = unlist(strsplit(cols[i,3], split = ""))[6:7]
      b = paste0(b[1],b[2])
      cols[i,4] = paste0("#",r,g,b)
    }
    
    m1 = cbind.data.frame(meta[, c("UMAP_1", "UMAP_2")], cols)
    #max1 = max(m1$UMAP_1)
    #min1 = min(m1$UMAP_1)
    
    #max2 = max(m1$UMAP_2)
    #min2 = min(m1$UMAP_2)
    
    #diff1 = max1 - min1
    #diff2 = max2 - min2
    
    
    par(bg="gray40")
    plot(m1$UMAP_1, m1$UMAP_2, pch = 20, cex = cex, cex.lab = 5, 
         ylab = "", main = paste("Level", nLev, m[1], m[2]), 
         col = m1$fcolor, cex.main = 5, 
         xlim = c(min1, max1 + 3))
    rasterImage(frink, max1, min2, max1 + (diff1/5), min2 + (diff2/5))
    
  } else {
    m = m[1:3]
    p1 = ggplot(meta, aes_string(x = "UMAP_1", y = "UMAP_2", color = m[1])) + 
      geom_point(size = 0.01)  + 
      theme(plot.title = element_text(size = 30)) + 
      scale_color_gradientn(colors = c("black", "white"))
    
    p2 = ggplot(meta, aes_string(x = "UMAP_1", y = "UMAP_2", color = m[2])) + 
      geom_point(size = 0.01) + 
      theme(plot.title = element_text(size = 30)) + 
      scale_color_gradientn(colors = c("black", "white"))
    
    p3 = ggplot(meta, aes_string(x = "UMAP_1", y = "UMAP_2", color = m[3])) + 
      geom_point(size = 0.01) + theme_classic() + 
      theme(plot.title = element_text(size = 30)) + 
      scale_color_gradientn(colors = c("black", "white"))
    
    cols = data.frame(layer_data(p1)$colour,
                      layer_data(p2)$colour,
                      layer_data(p3)$colour)
    cols$fcolor = NA
    for(i in 1:nrow(cols)){
      r = unlist(strsplit(cols[i,1], split = ""))[2:3]
      r = paste0(r[1],r[2])
      g = unlist(strsplit(cols[i,2], split = ""))[4:5]
      g = paste0(g[1],g[2])
      b = unlist(strsplit(cols[i,3], split = ""))[6:7]
      b = paste0(b[1],b[2])
      cols[i,4] = paste0("#",r,g,b)
    }
    m1 = cbind.data.frame(meta[, c("UMAP_1", "UMAP_2")], cols)
    #m1 = cbind.data.frame(meta[, c("UMAP_1", "UMAP_2")], cols)
    #max1 = max(m1$UMAP_1)
    #min1 = min(m1$UMAP_1)
    
    #max2 = max(m1$UMAP_2)
    #min2 = min(m1$UMAP_2)
    
    #diff1 = max1 - min1
    #diff2 = max2 - min2
    de
    par(bg="gray40")
    plot(m1$UMAP_1, m1$UMAP_2, pch = 20, cex = cex, cex.lab = 5, 
         ylab = "", xlab = fclus,
         main = paste("Level", nLev, m[1], m[2], m[3]), 
         col = m1$fcolor, cex.main = 5, 
         xlim = c(min1, max1 + 3))
    rasterImage(frink, max1, min2, max1 + (diff1/5), min2 + (diff2/5))
  }
}