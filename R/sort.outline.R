sort.outline <- function(outline) {
  outline <- as.data.frame(which(outline == 1, arr.ind=T))
  outline$order <- 0
  outline$order[1] <- 1
  nearest <- 1
  for (i in 1:nrow(outline)) {
    distance <- sqrt((outline[nearest,1] - outline[ ,1])^2 + (outline[nearest,2] - outline[ ,2])^2)
    distance[outline$order != 0] <- 100000000000
    nearest <- which(distance==min(distance))[1]
    outline$order[nearest] <- i
  }
  outline <- outline[order(outline$order),1:2]
  return(outline)
}
