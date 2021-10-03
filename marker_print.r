marker_print = function(m1, m2){
  par(mfrow = c(1,2))
  features = m1
  plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n', cex.main = 5,
       xaxt='n', yaxt='n', xlab='', ylab='', main = "Final markers")
  pos1 = 0.1
  text = 1
  for(pos in length(features):1){
    if(pos %% 52 == 0){
      text = text + 2
      pos1 = 0.1
    }
    text(text, pos1, features[pos], pos = 1, cex = 5, col = "white")
    pos1 = pos1 + 0.1  
  }
  
  features = m2
  plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n', cex.main = 5,
       xaxt='n', yaxt='n', xlab='', ylab='', main = "Level wise markers")
  pos1 = 0.1
  text = 1
  for(pos in length(features):1){
    if(pos %% 52 == 0){
      text = text + 2
      pos1 = 0.1
    }
    text(text, pos1, features[pos], pos = 1, cex = 5, col = "black")
    pos1 = pos1 + 0.1  
  }
}