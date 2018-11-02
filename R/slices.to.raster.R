slices.to.raster <- function(in.array, n.row, n.col, orientation) {
  #-------------------------------------------------------------------------------------
  # Copyright (C) 2017 Koscik, Timothy R. All Rights Reserved
  #-------------------------------------------------------------------------------------
  
  slice.count <- 0
  out.mx <- numeric()
  for (i in 1:n.row) {
    row.temp <- numeric()
    for (j in 1:n.col) {
      slice.count <- slice.count + 1
      row.temp <- switch(orientation,
        `coronal`=rbind(row.temp, in.array[ ,slice.count, ]),
        `axial`=rbind(row.temp, in.array[ , ,slice.count]),
        `sagittal`=rbind(row.temp, in.array[slice.count, , ]))
    }
    out.mx <- cbind(row.temp, out.mx)
  }
  out.mx <- melt(out.mx)
  return(out.mx)
}