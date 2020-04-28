draw.ortho <- function(coords,
                       bg.nii,
                       bg.vol,
                       bg.mask=NULL,
                       bg.mask.vol=NULL,
                       bg.color=c("#000000", "#ffffff"),
                       bg.range=c(0.01,0.99),
                       bg.alpha=1,
                       fg.nii=NULL,
                       fg.vol=NULL,
                       fg.minmax=c(NULL,NULL),
                       fg.label="Intensity",
                       fg.mask=NULL,
                       fg.mask.vol=NULL,
                       fg.color=c("#440154", "#3C4984", "#26828B", "#49B570", "#B5D940",
                                  "#FDD626", "#FBB330", "#F38A47", "#E2665F", "#CC4678"),
                       fg.alpha=1,
                       roi.nii=NULL,
                       roi.vol=NULL,
                       roi.color="#ff00ff",
                       save.plot=F, save.dir, file.name,
                       img.format="png", img.w=17.6, img.unit="cm", img.dpi=600) {
  
  hdr <- nii.hdr(bg.nii)
  if (hdr$qform_code %in% c(1,-1)) {
    mm.coords <- world.to.mni(matrix(coords,ncol=3),
                              tform.xyz = matrix(c(hdr$pixdim[2],0,0,hdr$qoffset_x,
                                                   0,hdr$pixdim[3],0,hdr$qoffset_y,
                                                   0,0,hdr$pixdim[4],hdr$qoffset_z),
                                                 ncol=4, byrow=T))
  } else {
    mm.coords <- world.to.mni(matrix(coords,ncol=3),
                              nii.file = bg.nii)
  }
  
  img.dims <- nii.dims(bg.nii)[1:3]
  rel.dims <- img.dims / max(img.dims)
  rel.dims <- rel.dims[c(2,1,1)]
  
  # Set plot theme ---------------------------------------------------------------
  theme.obj <- theme(plot.title = element_blank(),
                     legend.position="none",
                     legend.title = element_blank(),
                     legend.text = element_text(size=8, margin=margin(1,0,0,0,"null")),
                     axis.title=element_blank(),
                     axis.text=element_blank(),
                     axis.ticks=element_blank(),
                     plot.subtitle = element_text(size=10, margin = margin(0,0,0,0,"null")),
                     plot.background=element_blank(),
                     panel.background=element_blank(),
                     panel.grid=element_blank(),
                     panel.border=element_blank(),
                     panel.spacing.x=unit(c(0,0,0,0),"null"),
                     panel.spacing.y=unit(c(0,0,0,0),"null"),
                     plot.margin=margin(0,0,0,0, "null"),
                     legend.margin=margin(0,0,0,0, "null"),
                     panel.spacing = margin(0,0,0,0, "null"))
  
  # Prep Background --------------------------------------------------------------
  nii.bg <- read.nii.volume(bg.nii, bg.vol)
  if (is.null(bg.mask)) {
    mask.bg <- array(1,dim=dim(nii.bg))
  } else {
    mask.bg <- read.nii.volume(bg.mask, bg.mask.vol)
  }
  x.mask <- mask.bg[coords[1], , ]
  x.mask[x.mask==0] <- NA
  x.mask[!is.na(x.mask)] <- 1
  x.bg <- melt(nii.bg[coords[1], , ] * x.mask)
  y.mask <- mask.bg[ , coords[2], ]
  y.mask[y.mask==0] <- NA
  y.mask[!is.na(y.mask)] <- 1
  y.bg <- melt(nii.bg[ , coords[2], ] * y.mask)
  z.mask <- mask.bg[ , , coords[3]]
  z.mask[z.mask==0] <- NA
  z.mask[!is.na(z.mask)] <- 1
  z.bg <- melt(nii.bg[ , , coords[3]] * z.mask)
  bg.all.vals <-c(x.bg$value, y.bg$value, z.bg$value)
  min.bg <- quantile(bg.all.vals, bg.range[1], na.rm=T)
  max.bg <- quantile(bg.all.vals, bg.range[2], na.rm=T)
  
  x.bg$scaled <- ((x.bg$value - min.bg)/(max.bg - min.bg)) * (diff(bg.range))
  x.bg$scaled[x.bg$scaled < 0] <- 0
  x.bg$scaled[x.bg$scaled > 1] <- 1
  temp.color <- cbind(colorRamp(bg.color)(x.bg$scaled),NA)
  temp.color[!is.na(x.bg$scaled),4] <- 255 * bg.alpha
  temp.color[is.na(x.bg$scaled), ] <- 0
  x.bg$color=rgb(red = temp.color[ ,1],
                 green = temp.color[ ,2],
                 blue = temp.color[ ,3],
                 alpha = temp.color[ ,4],
                 maxColorValue = 255)
  
  y.bg$scaled <- (y.bg$value - min.bg)/(max.bg - min.bg)
  y.bg$scaled[y.bg$scaled < 0] <- 0
  y.bg$scaled[y.bg$scaled > 1] <- 1
  temp.color <- cbind(colorRamp(bg.color)(y.bg$scaled),NA)
  temp.color[!is.na(y.bg$scaled),4] <- 255 * bg.alpha
  temp.color[is.na(y.bg$scaled), ] <- 0
  y.bg$color=rgb(red = temp.color[ ,1],
                 green = temp.color[ ,2],
                 blue = temp.color[ ,3],
                 alpha = temp.color[ ,4],
                 maxColorValue = 255)
  
  z.bg$scaled <- (z.bg$value - min.bg)/(max.bg - min.bg)
  z.bg$scaled[z.bg$scaled < 0] <- 0
  z.bg$scaled[z.bg$scaled > 1] <- 1
  temp.color <- cbind(colorRamp(bg.color)(z.bg$scaled),NA)
  temp.color[!is.na(z.bg$scaled),4] <- 255 * bg.alpha
  temp.color[is.na(z.bg$scaled), ] <- 0
  z.bg$color=rgb(red = temp.color[ ,1],
                 green = temp.color[ ,2],
                 blue = temp.color[ ,3],
                 alpha = temp.color[ ,4],
                 maxColorValue = 255)
  
  # Plot background --------------------------------------------------------------
  plot.x <- ggplot() +
    theme_bw() +
    coord_equal(expand=FALSE, clip="off") +
    geom_raster(data=x.bg, aes(x=Var1, y=Var2, fill=value), fill=x.bg$color) +
    annotate("text", x=round(median(x.bg$Var1)), y=-Inf, label=round(mm.coords[1],digits=2), vjust=0, hjust=0.5) +
    theme.obj
  
  plot.y <- ggplot() +
    theme_bw() +
    coord_equal(expand=FALSE, clip="off") +
    geom_raster(data=y.bg, aes(x=Var1, y=Var2, fill=value), fill=y.bg$color) +
    annotate("text", x=round(median(y.bg$Var1)), y=-Inf, label=round(mm.coords[2],digits=2), vjust=0, hjust=0.5) +
    annotate("text", x=-Inf, y=Inf, label="R", vjust=0, hjust=0) +
    annotate("text", x=Inf, y=Inf, label="L", vjust=0, hjust=1) +
    theme.obj
  
  plot.z <- ggplot() +
    theme_bw() +
    coord_equal(expand=FALSE, clip="off") +
    geom_raster(data=z.bg, aes(x=Var1, y=Var2, fill=value), fill=z.bg$color) +
    annotate("text", x=round(median(z.bg$Var1)), y=-Inf, label=round(mm.coords[3],digits=2), vjust=0, hjust=0.5) +
    theme.obj
  
  # Prep Foreground --------------------------------------------------------------
  if (!is.null(fg.nii)) {
    nii.fg <- read.nii.volume(fg.nii, fg.vol)
    mask.fg <- read.nii.volume(fg.mask, fg.mask.vol)
    
    x.mask <- mask.fg[coords[1], , ]
    x.mask[x.mask==0] <- NA
    x.mask[!is.na(x.mask)] <- 1
    x.fg <- melt(nii.fg[coords[1], , ] * x.mask)
    
    y.mask <- mask.fg[ , coords[2], ]
    y.mask[y.mask==0] <- NA
    y.mask[!is.na(y.mask)] <- 1
    y.fg <- melt(nii.fg[ , coords[2], ] * y.mask)
    
    z.mask <- mask.fg[ , , coords[3]]
    z.mask[z.mask==0] <- NA
    z.mask[!is.na(z.mask)] <- 1
    z.fg <- melt(nii.fg[ , , coords[3]] * z.mask)
    all.fg.vals <- c(x.fg$value, y.fg$value, z.fg$value)
    
    if (is.null(fg.minmax[1])) {
      min.fg <- min(all.fg.vals, na.rm=T)
    } else {
      min.fg <- fg.minmax[1]
    }
    if (is.null(fg.minmax[2])) {
      max.fg <- max(all.fg.vals, na.rm=T)
    } else {
      max.fg <- fg.minmax[2]
    }
    
    x.fg$scaled <- (x.fg$value - min.fg)/(max.fg - min.fg)
    temp.color <- cbind(colorRamp(fg.color)(x.fg$scaled),NA)
    temp.color[!is.na(x.fg$scaled),4] <- 255 * fg.alpha
    temp.color[is.na(x.fg$scaled), ] <- 0
    temp.color[temp.color > 255] <- 255
    temp.color[temp.color < 0] <- 0
    temp.color[is.na(temp.color)] <- 0
    x.fg$color=rgb(red = temp.color[ ,1],
                   green = temp.color[ ,2],
                   blue = temp.color[ ,3],
                   alpha = temp.color[ ,4],
                   maxColorValue = 255)
    
    y.fg$scaled <- (y.fg$value - min.fg)/(max.fg - min.fg)
    temp.color <- cbind(colorRamp(fg.color)(y.fg$scaled),NA)
    temp.color[!is.na(y.fg$scaled),4] <- 255 * fg.alpha
    temp.color[is.na(y.fg$scaled), ] <- 0
    temp.color[temp.color > 255] <- 255
    temp.color[temp.color < 0] <- 0
    temp.color[is.na(temp.color)] <- 0
    y.fg$color=rgb(red = temp.color[ ,1],
                   green = temp.color[ ,2],
                   blue = temp.color[ ,3],
                   alpha = temp.color[ ,4],
                   maxColorValue = 255)
    
    z.fg$scaled <- (z.fg$value - min.fg)/(max.fg - min.fg)
    temp.color <- cbind(colorRamp(fg.color)(z.fg$scaled),NA)
    temp.color[!is.na(z.fg$scaled),4] <- 255 * fg.alpha
    temp.color[is.na(z.fg$scaled), ] <- 0
    temp.color[temp.color > 255] <- 255
    temp.color[temp.color < 0] <- 0
    temp.color[is.na(temp.color)] <- 0
    z.fg$color=rgb(red = temp.color[ ,1],
                   green = temp.color[ ,2],
                   blue = temp.color[ ,3],
                   alpha = temp.color[ ,4],
                   maxColorValue = 255)
    
    ## Prep Foreground Legend
    if (is.null(fg.minmax[1])) {
      min.val <- min(c(x.fg$value, y.fg$value, z.fg$value), na.rm=T)
    } else {
      min.val <- min.fg
    }
    if (is.null(fg.minmax[2])) {
      max.val <- max(c(x.fg$value, y.fg$value, z.fg$value), na.rm=T)
    } else {
      max.val <- max.fg
    }
    
    plot.x <- plot.x +
      geom_raster(data=x.fg, aes(x=Var1, y=Var2, fill=value), fill=x.fg$color)
    plot.y <- plot.y +
      geom_raster(data=y.fg, aes(x=Var1, y=Var2, fill=value), fill=y.fg$color)
    plot.z <- plot.z +
      geom_raster(data=z.fg, aes(x=Var1, y=Var2, fill=value), fill=z.fg$color)
    
    ## Prep Foreground Legend
    fg.legend <- melt(matrix(seq(min.val, max.val, length.out = 100), ncol=10, nrow=100))
    plot.legend <- ggplot(fg.legend, aes(x=Var2, y=Var1, fill=value)) +
      theme_bw() +
      coord_equal(expand=FALSE, clip="off") +
      scale_fill_gradientn(colors=fg.color) +
      geom_raster() +
      theme.obj
  }
  
  # Prep ROI outlines ------------------------------------------------------------
  if (!is.null(roi.nii)) {
    n.roi <- length(roi.nii)
    for (i in 1:n.roi) {
      x.idx <- melt(read.nii.volume(roi.nii[i],1)[coords[1], , ])
      x.idx <- x.idx[x.idx$value > 0, ]
      x.roi <- fortify(rasterToPolygons(rasterFromXYZ(x.idx), dissolve=T))
      plot.x <- plot.x +
        geom_path(data=x.roi, aes(x=long, y=lat, group=group),
                  size=0.25, color=roi.color[i])
      
      y.idx <- melt(read.nii.volume(roi.nii[i],1)[ , coords[2], ])
      y.idx <- y.idx[y.idx$value > 0, ]
      y.roi <- fortify(rasterToPolygons(rasterFromXYZ(y.idx), dissolve=T))
      plot.y <- plot.y +
        geom_path(data=y.roi, aes(x=long, y=lat, group=group),
                  size=0.25,  color=roi.color[i])
      
      z.idx <- melt(read.nii.volume(roi.nii[i],1)[ , , coords[3]])
      z.idx <- z.idx[z.idx$value > 0, ]
      z.roi <- fortify(rasterToPolygons(rasterFromXYZ(z.idx), dissolve=T))
      plot.z <- plot.z +
        geom_path(data=z.roi, aes(x=long, y=lat, group=group),
                  size=0.25,  color=roi.color[i])
    }
  }
  
  if (is.null(fg.nii)) {
    the.plot <- arrangeGrob(plot.x, plot.y, plot.z,
                            nrow=1, widths=c(rel.dims),
                            respect=TRUE, clip=T, padding=unit(0.5,"lines"),
                            right=fg.label)
  } else {
    the.plot <- arrangeGrob(plot.x, plot.y, plot.z,
                            rectGrob(gp=gpar(col="white")),
                            arrangeGrob(plot.legend, nrow=1,
                                        top=as.character(round(max.val,2)),
                                        bottom=as.character(round(min.val,2))),
                            nrow=1, widths=c(rel.dims,0.05, 0.1),
                            respect=TRUE, clip=T, padding=unit(0.5,"lines"),
                            right=fg.label)
  }

  if (save.plot) {
    ggsave(filename = paste(save.dir, file.name, sep="/"),
           plot = the.plot,
           device = img.format,
           width = img.w,
           height = img.w * (rel.dims[1]/sum(c(rel.dims, 0.1))),
           units = img.unit,
           dpi = img.dpi)
  }
  return(the.plot)
}
