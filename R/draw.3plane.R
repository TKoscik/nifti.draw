draw.3plane <- function(anat.nii,
                        coords,
                        slice.order=c(1,2,3),
                        slice.rot90 = c(0,1,0),
                        over.nii, over.vol, over.color,
                        mask.nii, mask.vol,
                        save.dir, file.name,
                        img.format="pdf", img.w=8.5, img.unit="cm", img.dpi=600,
                        label.size=1,
                        save.plot=TRUE, return.plot=FALSE) {

  rotate <- function(x) t(apply(x, 2, rev))
#
  n.check <- c(length(over.nii), length(over.color), length(mask.nii))
  n.overlays <- max(n.check)
  if (!all(n.check == n.overlays)) {
    if (n.check[1] != n.overlays) {
      temp <- over.nii
      over.nii <- vector("list", n.overlays)
      for (i in 1:n.overlays) { over.nii[[i]] <- temp }
      over.vol <- rep(over.vol, n.overlays)
    }
    if (n.check[2] != n.overlays) {
      temp <- over.color
      over.color <- vector("list", n.overlays)
      for (i in 1:n.overlays) { over.color[[i]] <- temp }
    }
    if (n.check[3] != n.overlays) {
      temp <- mask.nii
      mask.nii <- vector("list", n.overlays)
      for (i in 1:n.overlays) { mask.nii[[i]] <- temp }
      temp <- mask.vol
      mask.vol <- vector("list", n.overlays)
      for (i in 1:n.overlays) { mask.vol[[i]] <- temp }
    }
  }

  if (missing(img.format)) { img.format <- "png" }

  # Load Anatomical ----
  img.anat <- read.nii.volume(anat.nii, 1)
  img.anat[img.anat==0] <- NA
  img.dims <- dim(img.anat)

  # Load Overlays ----
  img.over <- vector("list", length=n.overlays)
  for (i in 1:n.overlays) { img.over[[i]] <- read.nii.volume(over.nii[[i]], over.vol[i]) }

  # Load masks ----
  img.mask <- vector("list", length=n.overlays)
  for (i in 1:n.overlays) {
    if (mask.nii[[i]] == "none") {
      img.mask[[i]] <- array(as.numeric(img.over[[i]] != 0), dim=dim(img.anat))
    } else {
      if (mask.vol[[i]] == "all") {
        mask.vols <- 1:(nii.dims(mask.nii[[i]])[4])
      } else if (is.numeric(mask.vol[[i]])) {
        mask.vols <- mask.vol[[i]]
      } else { stop("Cannot parse mask volumes") }

      img.mask[[i]] <- array(0, dim=dim(img.anat))
      for (j in mask.vols) {
        img.mask[[i]] <- img.mask[[i]] + read.nii.volume(mask.nii[[i]], j)
      }
      img.mask[[i]][img.mask[[i]]>1] <- 1
    }
    img.over[[i]][img.mask[[i]]==0] <- NA
  }

# Anatomical Slices ------------------------------------------------------------
  slice1 <- switch(as.character(slice.order[1]),
                   `1`=img.anat[coords[1], , ],
                   `2`=img.anat[ , coords[1], ],
                   `3`=img.anat[ , , coords[1]])
  slice2 <- switch(as.character(slice.order[2]),
                   `1`=img.anat[coords[2], , ],
                   `2`=img.anat[ , coords[2], ],
                   `3`=img.anat[ , , coords[2]])
  slice3 <- switch(as.character(slice.order[3]),
                   `1`=img.anat[coords[3], , ],
                   `2`=img.anat[ , coords[3], ],
                   `3`=img.anat[ , , coords[3]])

  if (dim(slice1)[1] != min(dim(slice1))) { slice1 <- rotate(slice1) }
  if (dim(slice2)[1] != dim(slice1)[1]) { slice2 <- rotate(slice2) }
  if (dim(slice3)[1] != dim(slice1)[1]) { slice3 <- rotate(slice3) }

  slice1 <- switch(as.character(slice.rot90[1]),
                   `0`=slice1, `1`=rotate(slice1), `2`=rotate(rotate(slice1)), `3`=rotate(rotate(rotate(slice1))))
  slice2 <- switch(as.character(slice.rot90[2]),
                   `0`=slice2, `1`=rotate(slice2), `2`=rotate(rotate(slice2)), `3`=rotate(rotate(rotate(slice2))))
  slice3 <- switch(as.character(slice.rot90[3]),
                   `0`=slice3, `1`=rotate(slice3), `2`=rotate(rotate(slice3)), `3`=rotate(rotate(rotate(slice3))))

  ## slice L R gets swapped somehow
  slice2 <- slice2[ , seq(ncol(slice2), 1, -1)]
  slice3 <- slice3[seq(nrow(slice3), 1, -1), ]

  anat.raster <- melt(cbind(slice1, slice2, slice3))

# MNI Labels --------------------------------------------------------------------
  mni.coords <- world.to.mni(matrix(coords,ncol=3), nii.file=anat.nii)
  mni.labels <- round(data.frame(xvar = c(dim(slice1)[2]/2,
                                    dim(slice1)[2] + dim(slice2)[2]/2,
                                    dim(slice1)[2] + dim(slice2)[2] + dim(slice3)[2]/2),
                           yvar = rep(dim(slice1)[1], 3),
                           labels = mni.coords), 1)

# set text height -------
  full.width <- dim(slice1)[2] + dim(slice2)[2] + dim(slice3)[2]
  full.height <- dim(slice1)[1] + n.overlays * 2 * round(dim(slice1)[1] * 0.05)
  wh.ratio <- full.width/full.height
  img.h <- img.w / wh.ratio
  label.pos <- round(dim(slice1)[1] * 0.05)
  label.height <- round(img.h) * 0.05 * label.size
  print(label.height)
  label.height <- switch(img.unit,
                         `cm`=label.height*10,
                         `in`=label.height*25.4)

# Overlay Slices ---------------------------------------------------------------
  over.raster <- vector("list", n.overlays)
  for (i in 1:n.overlays) {
    slice1 <- switch(as.character(slice.order[1]),
                     `1`=img.over[[i]][coords[1], , ],
                     `2`=img.over[[i]][ , coords[1], ],
                     `3`=img.over[[i]][ , , coords[1]])
    slice2 <- switch(as.character(slice.order[2]),
                     `1`=img.over[[i]][coords[2], , ],
                     `2`=img.over[[i]][ , coords[2], ],
                     `3`=img.over[[i]][ , , coords[2]])
    slice3 <- switch(as.character(slice.order[3]),
                     `1`=img.over[[i]][coords[3], , ],
                     `2`=img.over[[i]][ , coords[3], ],
                     `3`=img.over[[i]][ , , coords[3]])

    mask1 <- switch(as.character(slice.order[1]),
                     `1`=img.mask[[i]][coords[1], , ],
                     `2`=img.mask[[i]][ , coords[1], ],
                     `3`=img.mask[[i]][ , , coords[1]])
    mask2 <- switch(as.character(slice.order[2]),
                     `1`=img.mask[[i]][coords[2], , ],
                     `2`=img.mask[[i]][ , coords[2], ],
                     `3`=img.mask[[i]][ , , coords[2]])
    mask3 <- switch(as.character(slice.order[3]),
                     `1`=img.mask[[i]][coords[3], , ],
                     `2`=img.mask[[i]][ , coords[3], ],
                     `3`=img.mask[[i]][ , , coords[3]])

    slice1 <- slice1 * mask1
    slice2 <- slice2 * mask2
    slice3 <- slice3 * mask3

    if (dim(slice1)[1] != min(dim(slice1))) { slice1 <- rotate(slice1) }
    if (dim(slice2)[1] != dim(slice1)[1]) { slice2 <- rotate(slice2) }
    if (dim(slice3)[1] != dim(slice1)[1]) { slice3 <- rotate(slice3) }

    slice1 <- switch(as.character(slice.rot90[1]),
                     `0`=slice1, `1`=rotate(slice1), `2`=rotate(rotate(slice1)), `3`=rotate(rotate(rotate(slice1))))
    slice2 <- switch(as.character(slice.rot90[2]),
                     `0`=slice2, `1`=rotate(slice2), `2`=rotate(rotate(slice2)), `3`=rotate(rotate(rotate(slice2))))
    slice3 <- switch(as.character(slice.rot90[3]),
                     `0`=slice3, `1`=rotate(slice3), `2`=rotate(rotate(slice3)), `3`=rotate(rotate(rotate(slice3))))
    ## slice L R gets swapped somehow
    slice2 <- slice2[ , seq(ncol(slice2), 1, -1)]
    slice3 <- slice3[seq(nrow(slice3), 1, -1), ]

    over.raster[[i]] <- melt(cbind(slice1, slice2, slice3))
  }

# Convert to raster matrix -----------------------------------------------------
  for (i in 1:n.overlays) {
    col.ls <- colorRamp(over.color[[i]])

    over.raster[[i]]$value.scale <- (over.raster[[i]]$value - min(over.raster[[i]]$value, na.rm=TRUE)) /
      (max(over.raster[[i]]$value, na.rm=TRUE) - min(over.raster[[i]]$value, na.rm=TRUE))
    color.temp <- col.ls(over.raster[[i]]$value.scale[which(!is.na(over.raster[[i]]$value.scale))])/255

    over.raster[[i]]$red <- numeric(length=nrow(over.raster[[i]]))
    over.raster[[i]]$green <- numeric(length=nrow(over.raster[[i]]))
    over.raster[[i]]$blue <- numeric(length=nrow(over.raster[[i]]))
    over.raster[[i]]$red[which(!is.na(over.raster[[i]]$value.scale))] <- over.raster[[i]]$red[which(!is.na(over.raster[[i]]$value.scale))] + color.temp[, 1]
    over.raster[[i]]$green[which(!is.na(over.raster[[i]]$value.scale))] <- over.raster[[i]]$green[which(!is.na(over.raster[[i]]$value.scale))] + color.temp[, 2]
    over.raster[[i]]$blue[which(!is.na(over.raster[[i]]$value.scale))] <- over.raster[[i]]$blue[which(!is.na(over.raster[[i]]$value.scale))] + color.temp[, 3]

    over.raster[[i]]$red[is.na(over.raster[[i]]$red)] <- 0
    over.raster[[i]]$red[over.raster[[i]]$red > 1] <- 1
    over.raster[[i]]$blue[is.na(over.raster[[i]]$blue)] <- 0
    over.raster[[i]]$blue[over.raster[[i]]$blue > 1] <- 1
    over.raster[[i]]$green[is.na(over.raster[[i]]$green)] <- 0
    over.raster[[i]]$green[over.raster[[i]]$green > 1] <- 1
    over.raster[[i]]$alpha <- as.numeric(over.raster[[i]]$red > 0) + as.numeric(over.raster[[i]]$green > 0) + as.numeric(over.raster[[i]]$blue > 0)
    over.raster[[i]]$alpha <- as.numeric(over.raster[[i]]$alpha > 0)
  }

  plotf <- anat.raster
  plotf$r <- numeric(nrow(plotf))
  plotf$g <- numeric(nrow(plotf))
  plotf$b <- numeric(nrow(plotf))
  plotf$a <- numeric(nrow(plotf))
  for (i in 1:n.overlays) {
    plotf$r <- plotf$r + over.raster[[i]]$red
    plotf$g <- plotf$g + over.raster[[i]]$green
    plotf$b <- plotf$b + over.raster[[i]]$blue
    plotf$a <- plotf$a + over.raster[[i]]$alpha
  }
  plotf$a[plotf$a > 1] <- 1

  plot.img <- ggplot(plotf, aes(x=Var2, y=Var1, fill=value)) +
    geom_raster() +
    geom_raster(inherit.aes=FALSE,
                data=plotf, aes(x=Var2, y=Var1),
                fill=rgb(red=plotf$r,
                         green=plotf$g,
                         blue=plotf$b,
                         alpha=plotf$a*1),
                na.rm=TRUE) +
    coord_equal() +
    scale_fill_gradient(low="#000000", high="#ffffff", na.value="transparent") +
    theme(legend.position="none",
          legend.spacing=unit(0,"null"),
          axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          plot.margin=unit(c(0,0,0,0), "null"),
          plot.background=element_rect(fill = "transparent",colour = NA),
          panel.spacing=unit(c(0,0,0,0), "null"),
          panel.background=element_rect(fill = "transparent",colour = NA),
          panel.grid=element_blank(),
          panel.border=element_blank()) +
    geom_text(inherit.aes = FALSE,
              data = mni.labels, aes(x=xvar, y=yvar, label = labels),
              color = "#484848", size=label.height)

# color bar labels
  cbar.labels <- data.frame(x=c(rep(mni.labels$xvar[1], n.overlays),
                                rep(mni.labels$xvar[3], n.overlays)),
                            y=seq(-label.pos,n.overlays*(-(label.pos*2)), -(label.pos*2)),
                            labels=numeric(n.overlays))
  for (i in 1:n.overlays) {
    cbar.labels$labels[i] <- round(min(over.raster[[i]]$value, na.rm=TRUE), digits=3)
    cbar.labels$labels[i+n.overlays] <- round(max(over.raster[[i]]$value, na.rm=TRUE), digits=3)
  }
  plot.img <- plot.img +
    geom_text(inherit.aes = FALSE,
              data = cbar.labels, aes(x=x, y=y, label = labels),
              color = "#484848", size=label.height, hjust="outward")

# color bars ------------------------
  cbar.raster <- vector("list", n.overlays)
  for (i in 1:n.overlays) {
    cbar.fcn <- colorRampPalette(over.color[[i]])
    cbar.raster[[i]] <- matrix(cbar.fcn(50^2), nrow=50)
    plot.img <- plot.img +
      annotation_raster(cbar.raster[[i]],
                        mni.labels$xvar[1]+1, mni.labels$xvar[3]-1,
                        ((2*label.pos*(i-1))+(label.pos/2))*(-1),
                        ((2*label.pos*(i-1))+(3*label.pos/2))*(-1))
  }

  if (save.plot) {
    ggsave(filename=paste0(save.dir, "/", file.name, ".", img.format),
           plot.img,
           width=img.w, height=img.h, unit=img.unit, dpi=img.dpi)
  }
  if (return.plot) { return(plot.img) }

}
