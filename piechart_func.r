piechart_func = function(i1, sample_percent){
  sl = as.data.frame(table(i1$sample_location))
  sl[,3] = (sl[,2]/sum(sl[,2]))*100
  sl[,4] = as.character(sl[,1])
  sl[,4] = sl[,4] %>%
    gsub("Core_SNU27_Core", "lightblue", .) %>%
    gsub("Core_SNU33_Core", "pink", .) %>% 
    gsub("Core_SNU33b2_Core", "hotpink", .) %>%
    gsub("Core_SNU40_Core", "springgreen", .) %>%
    gsub("Core_SNU43_Core", "yellow", .) %>%
    gsub("Accutase_UK", "purple", .) %>%
    gsub("GBM_v7_UK", "mediumorchid4", .) %>%
    gsub("Peri_SNU27_Peri", "navy", .) %>%
    gsub("Peri_SNU33_Peri", "hotpink3", .) %>%
    gsub("Peri_SNU33b2_Peri", "hotpink4", .) %>%
    gsub("Peri_SNU40_Peri", "olivedrab", .) %>%
    gsub("Peri_SNU43_Peri", "goldenrod1", .) %>%
    gsub("SNU38_Control", "black", .)
  
  colnames(sl) = c("location", "val", "percent", "color")
  rownames(sl) = sl[,1]
  sl = merge(sl, sample_percent, all.x = T, by = 0)
  sl = sl[,-1]
  sl$spercent = (sl$val/sl$freq) * 100
  sl$label = paste0(sl$val, "(", round(sl$spercent, 1), "%)")
  
  num1 = grep("Core", sl[,1])
  num2 = grep("Peri", sl[,1])
  
  mat1 = sl[num1,]
  mat2 = sl[num2,]
  mat1 = mat1[, c("val", "freq")]
  mat2 = mat2[, c("val", "freq")]
  mat = matrix(NA, nrow = 2, ncol = 2)
  rownames(mat) = c("Core", "Peri")
  colnames(mat) = c("Present", "Absent")
  
  mat["Core","Present"] = sum(mat1[,1])
  mat["Core","Absent"] = sum(mat1[,2]) - sum(mat1[,1])
  
  mat["Peri","Present"] = sum(mat2[,1])
  mat["Peri","Absent"] = sum(mat2[,2]) - sum(mat2[,1])
  
  if(max(mat)>=10){
    a = chisq.test(mat)
  
    pie(sl$percent, sl$label, col = sl$color, cex = 5,
        main = paste0("Chisq pvalue: ", round(a$p.value, 3)), cex.main = 5)
    legend("topright", legend = sl$location, col = sl$color, 
           pch = 20, cex = 3)
    return(a$p.value)
  }else{
    pie(sl$percent, sl$label, col = sl$color, cex = 5,
        main = paste0("Chisq pvalue: ", "NA"), cex.main = 5)
    legend("topright", legend = sl$location, col = sl$color, 
           pch = 20, cex = 3)
    return (TRUE)
  }
  
}