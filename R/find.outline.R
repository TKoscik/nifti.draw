find.outline <- function(mask) {
  outline <- mask * 0
  idx <- which(mask != 0, arr.ind=T)
  # neighborhood <- matrix(c(-1,1,-1,0,-1,1,0,-1,0,1,1,-1,1,0,1,1),ncol=2, byrow=T)
  neighborhood <- matrix(c(-1,0,0,-1,0,1,1,0), ncol=2, byrow=T)
  for (i in 1:nrow(idx)) {
    new.idx <- cbind(idx[i,1] + neighborhood[ ,1],
                     idx[i,2] + neighborhood[ ,2])
    if (all(mask[new.idx] != 0)) {
      outline[idx[i,1],idx[i,2]] <- 0
    } else {
      outline[idx[i,1],idx[i,2]] <- 1
    }
  }
  return(outline)
}
