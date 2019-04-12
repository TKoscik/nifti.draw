draw.overlay <- function(anat.nii,
                         over.nii, over.vol, over.color,
                         mask.nii, mask.vol,
                         roi.nii=NULL, roi.val=NULL, roi.color="#ff64ff",
                         orientation = "coronal",
                         idx.slice=TRUE, cbars=TRUE,
                         save.dir, file.name,
                         img.format="png", img.w=NULL, img.unit="cm", img.dpi=600,
                         save.plot=FALSE, return.plot=TRUE) {

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

  # Load mask ----
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

  # Load ROIs ----
  if (!is.null(roi.nii)) {
  img.roi <- read.nii.volume(roi.nii, 1)
  ex.rois <- c(0, which(!(1:max(img.roi) %in% roi.val)))
  for (i in ex.rois) { img.roi[img.roi==i] <- NA }
  }

  # Get Slices
  slices <- numeric(0)
  if (!is.null(roi.nii)) {
    all.mask <- !is.na(img.roi)
  } else {
    all.mask <- array(0, dim=dim(img.over[[1]]))
  }
  for (i in 1:n.overlays) {
    all.mask <- all.mask + !is.na(img.over[[i]])
  }

  if (orientation=="coronal" | orientation=="c") {
    for (i in 1:dim(all.mask)[2]) {
      if (sum(all.mask[ ,i, ]) > 0) {
        slices <- c(slices, i)
      }
    }
    slices <- slices[-c(1,length(slices))]
    slice.length <- c(Inf,40,35,30,24,20,15,12,6,4,3,2,1)
    slice.length <- slice.length[min(which((length(slices) >= slice.length)==TRUE))]
    n.row <- switch(as.character(slice.length),
                    `40`=5, `35`=5, `30`=5, `24`=4, `20`=4, `15`=3, `12`=3, `6`=2,
                    `4`=2, `3`=3, `2`=2, `1`=1)
    n.col <- switch(as.character(slice.length),
                    `40`=8, `35`=7, `30`=6, `24`=6, `20`=5, `15`=5, `12`=4, `6`=3,
                    `4`=2, `3`=1, `2`=1, `1`=1)
    slices <- slices[round(seq(from=1, to=length(slices), length.out=slice.length))]

    img.idx <- img.anat[floor(dim(img.anat)[1]/2), , ]
    img.idx <- melt(img.idx)

    img.anat <- img.anat[ ,slices, ]
    img.anat <- slices.to.raster(img.anat,n.row,n.col,"coronal")
    for (i in 1:n.overlays) {
      img.over[[i]] <- img.over[[i]][ ,slices, ]
      img.over[[i]] <- slices.to.raster(img.over[[i]],n.row,n.col,"coronal")
    }
    if (!is.null(roi.nii)) {
      img.roi <- img.roi[ ,slices, ]
      img.roi <- slices.to.raster(img.roi,n.row,n.col,"coronal")
    }


  } else if (orientation=="axial" | orientation=="a") {
    for (i in 1:dim(all.mask)[3]) {
      if (sum(all.mask[ , ,i]) > 0) {
        slices <- c(slices, i)
      }
    }
    slices <- slices[-c(1,length(slices))]
    slice.length <- c(Inf,40,35,30,24,20,15,12,6,4,3,2,1)
    slice.length <- slice.length[min(which((length(slices) >= slice.length)==TRUE))]
    n.row <- switch(as.character(slice.length),
                    `40`=5, `35`=5, `30`=5, `24`=4, `20`=4, `15`=3, `12`=3, `6`=2,
                    `4`=2, `3`=3, `2`=2, `1`=1)
    n.col <- switch(as.character(slice.length),
                    `40`=8, `35`=7, `30`=6, `24`=6, `20`=5, `15`=5, `12`=4, `6`=3,
                    `4`=2, `3`=1, `2`=1, `1`=1)
    slices <- slices[round(seq(from=1, to=length(slices), length.out=slice.length))]

    img.idx <- img.anat[ ,floor(dim(img.anat)[2]/2), ]
    img.idx <- melt(img.idx)

    img.anat <- img.anat[ , ,slices]
    img.anat <- slices.to.raster(img.anat,n.row,n.col,"axial")
    for (i in 1:n.overlays) {
      img.over[[i]] <- img.over[[i]][ , ,slices]
      img.over[[i]] <- slices.to.raster(img.over[[i]],n.row,n.col,"axial")
    }
    if (!is.null(roi.nii)) {
      img.roi <- img.roi[ , ,slices]
      img.roi <- slices.to.raster(img.roi,n.row,n.col,"axial")
    }


  } else if (orientation=="sagittal" | orientation=="s") {
    for (i in 1:dim(all.mask)[1]) {
      if (sum(all.mask[i, , ]) > 0) {
        slices <- c(slices, i)
      }
    }
    slices <- slices[-c(1,length(slices))]
    slice.length <- c(Inf,40,35,30,24,20,15,12,6,4,3,2,1)
    slice.length <- slice.length[min(which((length(slices) >= slice.length)==TRUE))]
    n.row <- switch(as.character(slice.length),
                    `40`=5, `35`=5, `30`=5, `24`=4, `20`=4, `15`=3, `12`=3, `6`=2,
                    `4`=2, `3`=3, `2`=2, `1`=1)
    n.col <- switch(as.character(slice.length),
                    `40`=8, `35`=7, `30`=6, `24`=6, `20`=5, `15`=5, `12`=4, `6`=3,
                    `4`=2, `3`=1, `2`=1, `1`=1)
    slices <- slices[round(seq(from=1, to=length(slices), length.out=slice.length))]

    img.idx <- img.anat[ , ,floor(dim(img.anat)[3]/2)]
    img.idx <- melt(img.idx)

    img.anat <- img.anat[slices, , ]
    img.anat <- slices.to.raster(img.anat,n.row,n.col,"sagittal")
    for (i in 1:n.overlays) {
      img.over[[i]] <- img.over[[i]][slices, , ]
      img.over[[i]] <- slices.to.raster(img.over[[i]],n.row,n.col,"sagittal")
    }
    if (!is.null(roi.nii)) {
      img.roi <- img.roi[slices, , ]
      img.roi <- slices.to.raster(img.roi,n.row,n.col,"sagittal")
    }


  } else { stop("Cannot parse orientation.") }

  # get slice MNI coordinates ----
  mni.labels <- data.frame(yvar = sort(rep(seq(0,1-1/n.row, length.out = n.row), n.col), decreasing = TRUE),
                           xvar = rep(seq(0,1-1/n.col, length.out = n.col), n.row),
                           labels = numeric(n.row * n.col))
  if (orientation=="coronal" | orientation=="c") {
    tform <- unlist(nii.hdr(anat.nii, "srow_y"))
    mni.labels$xvar <- mni.labels$xvar * nii.dims(anat.nii)[1] * n.col + nii.dims(anat.nii)[1] / 2
    mni.labels$yvar <- mni.labels$yvar * nii.dims(anat.nii)[3] * n.row - 3
    mni.labels$labels <- (slices - 1) * tform[2] + tform[4]
  } else if (orientation == "axial" | orientation == "a") {
    tform <- unlist(nii.hdr(anat.nii, "srow_z"))
    mni.labels$xvar <- mni.labels$xvar * nii.dims(anat.nii)[1] * n.col + nii.dims(anat.nii)[1] / 2
    mni.labels$yvar <- mni.labels$yvar * nii.dims(anat.nii)[2] * n.row - 3
    mni.labels$labels <- (slices - 1) * tform[3] + tform[4]
  } else if (orientation == "sagittal" | orientation == "s") {
    tform <- unlist(nii.hdr(anat.nii, "srow_x"))
    mni.labels$xvar <- mni.labels$xvar * nii.dims(anat.nii)[2] * n.col + nii.dims(anat.nii)[2] / 2
    mni.labels$yvar <- mni.labels$yvar * nii.dims(anat.nii)[3] * n.row - 3
    mni.labels$labels <- (slices - 1) * tform[1] + tform[4]
  } else { stop("Cannot parse orientation.") }

  # Convert to raster matrix ----
  for (i in 1:n.overlays) {
    col.ls <- colorRamp(over.color[[i]])

    img.over[[i]]$value.scale <- (img.over[[i]]$value - min(img.over[[i]]$value, na.rm=TRUE)) /
      (max(img.over[[i]]$value, na.rm=TRUE) - min(img.over[[i]]$value, na.rm=TRUE))
    color.temp <- col.ls(img.over[[i]]$value.scale[which(!is.na(img.over[[i]]$value.scale))])/255

    img.over[[i]]$red <- numeric(length=nrow(img.over[[i]]))
    img.over[[i]]$green <- numeric(length=nrow(img.over[[i]]))
    img.over[[i]]$blue <- numeric(length=nrow(img.over[[i]]))
    img.over[[i]]$red[which(!is.na(img.over[[i]]$value.scale))] <- img.over[[i]]$red[which(!is.na(img.over[[i]]$value.scale))] + color.temp[, 1]
    img.over[[i]]$green[which(!is.na(img.over[[i]]$value.scale))] <- img.over[[i]]$green[which(!is.na(img.over[[i]]$value.scale))] + color.temp[, 2]
    img.over[[i]]$blue[which(!is.na(img.over[[i]]$value.scale))] <- img.over[[i]]$blue[which(!is.na(img.over[[i]]$value.scale))] + color.temp[, 3]

    img.over[[i]]$red[is.na(img.over[[i]]$red)] <- 0
    img.over[[i]]$red[img.over[[i]]$red > 1] <- 1
    img.over[[i]]$blue[is.na(img.over[[i]]$blue)] <- 0
    img.over[[i]]$blue[img.over[[i]]$blue > 1] <- 1
    img.over[[i]]$green[is.na(img.over[[i]]$green)] <- 0
    img.over[[i]]$green[img.over[[i]]$green > 1] <- 1
    img.over[[i]]$alpha <- as.numeric(img.over[[i]]$red > 0) + as.numeric(img.over[[i]]$green > 0) + as.numeric(img.over[[i]]$blue > 0)
    img.over[[i]]$alpha <- as.numeric(img.over[[i]]$alpha > 0)
  }

  plotf <- img.anat
  plotf$r <- numeric(nrow(plotf))
  plotf$g <- numeric(nrow(plotf))
  plotf$b <- numeric(nrow(plotf))
  plotf$a <- numeric(nrow(plotf))
  for (i in 1:n.overlays) {
    plotf$r <- plotf$r + img.over[[i]]$red
    plotf$g <- plotf$g + img.over[[i]]$green
    plotf$b <- plotf$b + img.over[[i]]$blue
    plotf$a <- plotf$a + img.over[[i]]$alpha
  }
  plotf$r[plotf$r > 1] <- 1
  plotf$g[plotf$g > 1] <- 1
  plotf$b[plotf$b > 1] <- 1
  plotf$a[plotf$a > 1] <- 1

  if (!is.null(roi.nii)) {
    roi.raster <- rasterFromXYZ(img.roi)
    roi.poly <- fortify(rasterToPolygons(roi.raster, dissolve=TRUE))
  }

  plot.img <- ggplot(plotf, aes(x=Var1, y=Var2, fill=value)) +
    geom_raster() +
    geom_raster(inherit.aes=FALSE,
                data=plotf, aes(x=Var1, y=Var2),
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
              color = "#484848", size=3)

  if (!is.null(roi.nii)) {
    plot.img <- plot.img + geom_path(inherit.aes=FALSE, data=roi.poly,
      aes(x=long, y=lat, group=group), size=0.25, alpha=0.5,
      color=roi.color, linetype="solid")
  }

  plot.idx <- ggplot(img.idx, aes(x=Var1, y=Var2, fill=value)) +
    geom_raster() +
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
          panel.border=element_blank(),
          text=element_text(size=14)) +
    annotate("text", x=-1, y=-2, label=1, size=n.row*.48)
  if (orientation == "coronal") {
    plot.idx <- plot.idx +
      annotate("segment", x=slices, xend=slices, y=0, yend=Inf, color="#0000ff", size=0.25) +
      annotate("text", x=img.dims[2]-4, y=-2, label=img.dims[2], size=n.row*.48) +
      annotate("segment", x=5, xend=img.dims[2]-15, y=-2, yend=-2, arrow=arrow(length=unit(.05,"npc")))
  } else if (orientation == "axial") {
    plot.idx <- plot.idx +
      annotate("segment", x=0, xend=Inf, y=slices, yend=slices, color="#0000ff", size=0.25) +
      annotate("text", y=img.dims[1]-4, x=-2, label=img.dims[1], size=n.row*.48) +
      annotate("segment", y=5, yend=img.dims[1]-15, x=-2, xend=-2, arrow=arrow(length=unit(.05,"npc")))
  } else if (orientation == "sagittal") {
    plot.idx <- plot.idx +
      annotate("segment", x=slices, xend=slices, y=0, yend=Inf, color="#0000ff", size=0.25) +
      annotate("text", x=img.dims[3]-4, y=-2, label=img.dims[3], size=n.row*.48) +
      annotate("segment", x=5, xend=img.dims[3]-15, y=-2, yend=-2, arrow=arrow(length=unit(.05,"npc")))
  }

  plot.cbar <- vector("list", n.overlays)
  scale.res <- 500
  img.cbar <- data.frame(x=1,y=1:scale.res)
  for (i in 1:n.overlays) {
    cbar.labels <- data.frame(xval=c(1,1),
                              yval=c(scale.res+0.2*scale.res, -0.2*scale.res),
                              the.labels=c(round(max(img.over[[i]]$value, na.rm=TRUE), digits=3),
                                           round(min(img.over[[i]]$value, na.rm=TRUE), digits=3)))
    plot.cbar[[i]] <- ggplot(img.cbar, aes(x=x, y=y, fill=y)) +
      geom_raster() +
      coord_equal(ratio=((n.col+1)/8)*0.05, expand=FALSE) +
      scale_fill_gradientn(colors=over.color[[i]]) +
      geom_text(inherit.aes=FALSE, data=cbar.labels,
                aes(x=xval, y=yval, label=the.labels),
                angle=90, size=n.row*.5) +
      ylim(c(-0.4*scale.res, scale.res+0.4*scale.res)) +
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
            panel.border=element_blank(),
            text=element_text(size=14))
  }

  lay <- cbind(matrix(1, nrow=n.row*n.overlays, ncol=n.col*n.overlays),
               rbind(matrix(sort(rep(2:(n.overlays+1), n.overlays*(n.row-1))),
                            nrow=n.overlays*(n.row-1), ncol=n.overlays),
                     matrix(n.overlays+2, nrow=n.overlays, ncol=n.overlays)))

  plot.ls <- paste(c("plot.img, ", sprintf("plot.cbar[[%0.0f]], ",1:n.overlays), "plot.idx"),collapse=" ")
  plot.all <- eval(parse(text=sprintf("arrangeGrob(%s, layout_matrix=lay)",plot.ls)))
  # plot.all <- arrangeGrob(plot.img, plot.cbar[[1]], plot.cbar[[2]], plot.idx, layout_matrix=lay)

  # ggsave(filename="/home/tim/Documents/duplex.20160104/test.roi.overlay.coronal.png", plot.all,
  #        width=n.col+1, height=n.row, unit="in", dpi=300)

  if (is.null(img.w)) {
    img.w <- (n.col+1)*2
    img.h <- (n.row)*2
  } else {
    img.h <- (n.row * img.w)/(n.col + 1)
  }

  if (save.plot) {
    ggsave(filename=paste0(save.dir, "/", file.name, ".", img.format),
           plot.all,
           width=img.w, height=img.h, unit=img.unit, dpi=img.dpi)
  }
  if (return.plot) { return(plot.all) }

}
