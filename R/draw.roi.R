draw.roi <- function(anat.img, anat.vol=1, anat.color=c("#ffffff", "#000000"),
                     roi.img, roi.vol=1, roi.color=c("#a80000", "#0000e8", "#006400"),
                     orientation=c("sagittal", "axial", "coronal"),
                     n.row = 5, n.col = 4) {

  # Debug ----
  # anat.img <- "/Shared/koscikt_scratch/nifti_qc_scratch/sub-231_ses-328zk16wb6_site-00201_acq-sagMPRAGEPROMO_T1w_prep-denoise.nii"
  # anat.vol <- 1
  # anat.color <- c("#ffffff", "#000000")
  # roi.img <- c("/Shared/koscikt_scratch/nifti_qc_scratch/sub-231_ses-328zk16wb6_site-00201_prep-bex0ANTS.nii",
  #              "/Shared/koscikt_scratch/nifti_qc_scratch/sub-231_ses-328zk16wb6_site-00201_prep-bex0MALF.nii")
  # roi.vol <- c(1,1)
  # roi.color <- c("#a80000", "#0000e8", "#006400", "#a86400", "#0064e8", "#a800e8")
  # orientation <- "sagittal"
  # ----

  # set slice dimension ----
  which.dim <- switch(orientation[1], `axial`=3, `coronal`=2, `sagittal`=1)

  # load images ----
  base.img <- read.nii.volume(anat.img, anat.vol)

  n.roi <- length(roi.img)
  roi <- vector("list", n.roi)
  roi.all <- array(0, dim=dim(base.img))
  for (i in 1:n.roi) {
    roi[[i]] <- read.nii.volume(roi.img[i], roi.vol[i])
    roi.all <- ((roi.all + roi[[i]]) > 0)*1
  }

  # determine slices to plot ----
  which.slices <- numeric()
  for (i in 1:dim(base.img)[which.dim]) {
    slice = switch(orientation[1],
                   `axial`=roi.all[ , ,i],
                   `coronal`=roi.all[ ,i, ],
                   `sagittal`=roi.all[i, , ])
    if (sum(slice) != 0) {
      which.slices <- c(which.slices, i)
    }
  }

  n.slices <- prod(n.row, n.col)
  if (length(which.slices) < n.slices) {
    calc.row <- floor(n.slices / n.col)
    if (calc.row < 1) {
      n.row <- 1
      n.col <- n.slices
    } else {
      n.row <- calc.row
      n.slices <- prod(n.row, n.col)
    }
  }
  slices <- which.slices[seq(0,length(which.slices)+1,length.out=n.slices+2)[-c(1,n.slices+2)]]

  # resample images to desired slices ----
  if (orientation[1] == "axial") {
    base.img <- base.img[ , , slices]
    for (i in 1:n.roi) { roi[[i]] <- roi[[i]][ , , slices] }
  } else if (orientation[1] == "coronal") {
    base.img <- base.img[ , slices, ]
    for (i in 1:n.roi) { roi[[i]] <- roi[[i]][ , slices, ] }
  } else if (orientation[1] == "sagittal") {
    base.img <- base.img[slices, , ]
    for (i in 1:n.roi) { roi[[i]] <- roi[[i]][slices, , ] }
  }

  # images to rasters & polygons ----
  base.img <- slices.to.raster(base.img, n.row, n.col, orientation[1])
  for (i in 1:n.roi) {
    roi[[i]] <- slices.to.raster(roi[[i]], n.row, n.col, orientation[1])
    roi[[i]] <- rasterFromXYZ(roi[[i]])
    roi[[i]] <- fortify(rasterToPolygons(roi[[i]], fun=function(x){x != 0}, na.rm=TRUE, digits=1, dissolve=TRUE))
  }

  # create plot
  roi.plot <- ggplot(base.img, aes(x=Var1, y=Var2, fill=value)) +
    geom_raster() +
    coord_equal() +
    scale_fill_gradient(low=anat.color[1], high=anat.color[2], na.value="transparent") +
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
          panel.border=element_blank())
  for (i in 1:n.roi) {
    if (i > length(roi.color)) {
      roi.color <- c(roi.color, paste0("#", paste0(as.hexmode(sample(0:255,size=3,replace=TRUE)), collapse=""), collapse=""))
    }
    roi.plot <- roi.plot +
      geom_path(inherit.aes=FALSE,
                data=roi[[i]], aes(x=long, y=lat, group=group),
                size=0.25, alpha=0.5, color=roi.color[i], linetype="solid")
  }

  return(roi.plot)
}
