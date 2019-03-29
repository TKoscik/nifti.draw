draw.effect <- function(model,
                        effect.name,
                        axes=NULL, labels=NULL,
                        plot.colors=NULL,
                        save.dir, file.name,
                        img.format="png", img.w=10, img.unit="cm", img.dpi=600,
                        save.plot=TRUE, return.plot=FALSE) {
  
  if (is.null(plot.colors)) {
    plot.colors <- colorRampPalette(c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF", "#FDC328FF", "#F89441FF", "#E56B5DFF", "#CC4678FF"))
  }
  
  if (class(model) == "merModLmerTest") {
    tempf <<- model@frame
    dv <- unlist(strsplit(as.character(model@call), split = "[ ~]"))[2]
    tempf[paste0(dv, ".dummy")] <- tempf[dv]
    temp.ranef <- ranef(model)
    for (i in 1:length(temp.ranef)) {
      which.vars <- names(temp.ranef)[[i]]
      which.vars <- unlist(strsplit(which.vars, split = ":"))
      tempf[paste0("dummy", i)] <- do.call(paste, c(tempf[which.vars], sep=":"))
      tempf[paste0(dv, ".dummy")] <- tempf[paste0(dv, ".dummy")] -
        unlist(temp.ranef[[i]])[as.numeric(factor(unlist(tempf[paste0("dummy", i)]), levels=unique(unlist(tempf[paste0("dummy", i)]))))]
    }
    tempf[dv] <- tempf[paste0(dv, ".dummy")]
    model <- update(model, . ~ ., data=tempf)
  } else {
    tempf <<- model@frame
  }

  # determine what to plot where
  if (is.null(axes)) {
    ef.vars <- unlist(strsplit(effect.name, split=":"))
    ef.num <- length(ef.vars)
    var.type <- sapply(tempf[ef.vars], class)
    axes <- switch(as.character(ef.num),
                   `1` = c(ef.vars[1], "fit"),
                   `2` = c(ef.vars[1], "fit", ef.vars[2]),
                   `3` = c(ef.vars[1], "fit", ef.vars[2], ef.vars[3]),
                   `4` = c(ef.vars[1], "fit", ef.vars[2], ef.vars[3], ef.vars[4]),
                   stop("Cannot plot interactions larger than 4-way"))
  } else {
    ef.vars <- axes[-which(axes=="fit")]
    ef.num <- length(ef.vars)
    var.type <- sapply(tempf[ef.vars], class)
  }

  # set labels
  if (is.null(labels)) {
    s.temp <- unlist(strsplit(axes, " "))
    labels <- paste(toupper(substring(s.temp,1,1)), substring(s.temp,2), sep="")
  }

  # Get effect levels if not supplied by user ----
  effect.levels <- numeric(ef.num)
    for (i in 1:ef.num) {
      if (i == 1) {
        if (var.type[i]=="factor") {
          effect.levels[i] <- eval(parse(text=sprintf("nlevels(tempf$%s)", ef.vars[i])))
        } else if (var.type[i]=="numeric" | var.type[i]=="integer") {
          effect.levels[i] <- 100
        } else { stop("Cannot parse variable class") }
      } else {
        if (var.type[i]=="factor") {
          effect.levels[i] <- eval(parse(text=sprintf("nlevels(tempf$%s)", ef.vars[i])))
        } else if (var.type[i]=="numeric" | var.type[i]=="integer") {
          if (length(unlist(unique(tempf[ef.vars[i]]))) >= 3) {
            effect.levels[i] <- 3
          } else {
            effect.levels[i] <- length(unlist(unique(tempf[ef.vars[i]])))
          }
        } else { stop("Cannot parse variable class") }
      }
  }

  ef.str <- switch(as.character(ef.num),
    `1`=sprintf("%s%s=%0.0f%s",
                "as.data.frame(effect(effect.name, model, xlevels=list(",
                ef.vars[1], effect.levels[1], ")))"),
    `2`=sprintf("%s%s=%0.0f, %s=%0.0f)))",
                "as.data.frame(effect(effect.name, model, xlevels=list(",
               ef.vars[1], effect.levels[1],
               ef.vars[2], effect.levels[2], ")))"),
    `3`=sprintf("%s%s=%0.0f, %s=%0.0f, %s=%0.0f%s",
                "as.data.frame(effect(effect.name, model, xlevels=list(",
               ef.vars[1], effect.levels[1],
               ef.vars[2], effect.levels[2],
               ef.vars[3], effect.levels[3], ")))"),
    `4`=sprintf("%s%s=%0.0f, %s=%0.0f, %s=%0.0f, %s=%0.0f%s",
                "as.data.frame(effect(effect.name, model, xlevels=list(",
               ef.vars[1], effect.levels[1],
               ef.vars[2], effect.levels[2],
               ef.vars[3], effect.levels[3],
               ef.vars[4], effect.levels[4], ")))"))
  ef <- eval(parse(text=ef.str))

  if (length(axes) > 2) {
    if (var.type[axes[3]]=="numeric") {
      ef[axes[3]] <- signif(ef[axes[3]], digits=3)
    }
  }
  if (ef.num == 3) {
    which.column <- which(colnames(ef)==axes[4])
    colnames(ef)[which.column] <- labels[4]
  }
  if (ef.num == 4) {
    which.column <- which(colnames(ef)==axes[4])
    colnames(ef)[which.column] <- labels[4]
    which.column <- which(colnames(ef)==axes[5])
    colnames(ef)[which.column] <- labels[5]
  }
  
  plot.colors <- plot.colors(length(unlist(unique(ef[axes[3]]))))
                             
  if (var.type[1]=="factor") {
    plot.str <- switch(as.character(ef.num),
      `1`=sprintf("ggplot(ef, aes(x=%s, y=%s, ymin=lower, ymax=upper)) + theme_bw() + geom_pointrange(size=1, shape=18) + xlab('%s') + ylab('%s') + theme(axis.title=element_text(size=14), axis.text=element_text(size=12), axis.text.x=element_text(angle=90, hjust=0, vjust=0.5), legend.position='bottom', legend.title=element_text(size=14), legend.text=element_text(size=12))", axes[1], axes[2], labels[1], labels[2]),
      `2`=sprintf("ggplot(ef, aes(x=%s, y=%s, ymin=lower, ymax=upper, fill=factor(%s), color=factor(%s))) + theme_bw() + geom_pointrange(size=1, shape=18, position=position_dodge(width=0.5)) + scale_fill_manual(values=plot.colors, name='%s') + scale_color_manual(values=plot.colors, name='%s') + xlab('%s') + ylab('%s') + theme(axis.title=element_text(size=14), axis.text=element_text(size=12), axis.text.x=element_text(angle=90, hjust=0, vjust=0.5), legend.position='bottom', legend.title=element_text(size=14), legend.text=element_text(size=12))", axes[1], axes[2], axes[3], axes[3], labels[3], labels[3], labels[1], labels[2]),
      `3`=sprintf("ggplot(ef, aes(x=%s, y=%s, ymin=lower, ymax=upper, fill=factor(%s), color=factor(%s))) + theme_bw() + geom_pointrange(size=1, shape=18, position=position_dodge(width=0.5)) + scale_fill_manual(values=plot.colors, name='%s') + scale_color_manual(values=plot.colors, name='%s') + facet_grid(.~`%s`, labeller=label_both) + xlab('%s') + ylab('%s') + theme(title=element_text(size=12), axis.title=element_text(size=14), axis.text=element_text(size=12), axis.text.x=element_text(angle=90, hjust=0, vjust=0.5), strip.background=element_rect(fill='#ffffff', size=0.5), strip.text=element_text(size=12), legend.position='bottom', legend.title=element_text(size=14), legend.text=element_text(size=12))", axes[1], axes[2], labels[3], axes[3], axes[3], labels[3], labels[4], labels[1], labels[2]),
      `4`=sprintf("ggplot(ef, aes(x=%s, y=%s, ymin=lower, ymax=upper, fill=factor(%s), color=factor(%s))) + theme_bw() + geom_pointrange(size=1, shape=18, position=position_dodge(width=0.5)) + scale_fill_manual(values=plot.colors, name='%s') + scale_color_manual(values=plot.colors, name='%s') + facet_grid(`%s`~`%s`, labeller=label_both) + xlab('%s') + ylab('%s') + theme(title=element_text(size=12), axis.title=element_text(size=14), axis.text=element_text(size=12), axis.text.x=element_text(angle=90, hjust=0, vjust=0.5), strip.background=element_rect(fill='#ffffff', size=0.5), strip.text=element_text(size=12), legend.position='bottom', legend.title=element_text(size=14), legend.text=element_text(size=12))", axes[1], axes[2], axes[3], axes[3], labels[3], labels[3], labels[5], labels[4], labels[1], labels[2]))
  } else if (var.type[1]=="numeric" | var.type[1]=="integer") {
    plot.str <- switch(as.character(ef.num),
      `1`=sprintf("ggplot(ef, aes(x=%s, y=%s, ymin=lower, ymax=upper)) + theme_bw() + geom_ribbon(alpha=0.25, color='transparent') + geom_line(size=1) + xlab('%s') + ylab('%s') + theme(axis.title=element_text(size=14), axis.text=element_text(size=12), axis.text.x=element_text(angle=90, hjust=0, vjust=0.5), legend.position='bottom', legend.title=element_text(size=14), legend.text=element_text(size=12))", axes[1], axes[2], labels[1], labels[2]),
      `2`=sprintf("ggplot(ef, aes(x=%s, y=%s, ymin=lower, ymax=upper, fill=factor(%s), color=factor(%s))) + theme_bw() + scale_fill_manual(values=plot.colors, name='%s') + scale_color_manual(values=plot.colors, name='%s') + geom_ribbon(color='transparent', alpha=0.25) + geom_line(size=1) + xlab('%s') + ylab('%s') + theme(axis.title=element_text(size=14), axis.text=element_text(size=12), axis.text.x=element_text(angle=90, hjust=0, vjust=0.5), legend.position='bottom', legend.title=element_text(size=14), legend.text=element_text(size=12))", axes[1], axes[2], axes[3], axes[3], labels[3], labels[3], labels[1], labels[2]),
      `3`=sprintf("ggplot(ef, aes(x=%s, y=%s, ymin=lower, ymax=upper, fill=factor(%s), color=factor(%s))) + theme_bw() + scale_fill_manual(values=plot.colors, name='%s') + scale_color_manual(values=plot.colors, name='%s') + geom_ribbon(color='transparent', alpha=0.25) + geom_line(size=1) + facet_grid(.~`%s`, labeller=label_both) + xlab('%s') + ylab('%s') + theme(axis.title=element_text(size=14), axis.text=element_text(size=12), axis.text.x=element_text(angle=90, hjust=0, vjust=0.5), strip.background=element_rect(fill='#ffffff', size=0.5), strip.text=element_text(size=12), legend.position='bottom', legend.title=element_text(size=14), legend.text=element_text(size=12))", axes[1], axes[2], axes[3], axes[3], labels[3], labels[3], labels[4], labels[1], labels[2]),
      `4`=sprintf("ggplot(ef, aes(x=%s, y=%s, ymin=lower, ymax=upper, fill=factor(%s), color=factor(%s))) + theme_bw() + scale_fill_manual(values=plot.colors, name='%s') + scale_color_manual(values=plot.colors, name='%s') + geom_ribbon(color='transparent', alpha=0.25) + geom_line(size=1) + facet_grid(`%s`~`%s`, labeller=label_both) + xlab('%s') + ylab('%s') + theme(axis.title=element_text(size=14), axis.text=element_text(size=12), axis.text.x=element_text(angle=90, hjust=0, vjust=0.5), strip.background=element_rect(fill='#ffffff', size=0.5), strip.text=element_text(size=12), legend.position='bottom', legend.title=element_text(size=14), legend.text=element_text(size=12))", axes[1], axes[2], axes[3], axes[3], labels[3], labels[3], labels[5], labels[4], labels[1], labels[2]))
  }
  the.plot <- eval(parse(text=plot.str))

  if (save.plot) {
    img.h <- switch(as.character(ef.num),
      `1`=(img.w/6)*5,
      `2`=(img.w/6)*6,
      `3`=(img.w/(3*effect.levels[3]))*(6),
      `4`=(img.w/(3*effect.levels[4]))*(6*effect.levels[3]))
    ggsave(filename=paste0(save.dir, "/", file.name, ".", img.format),
           plot=the.plot, width=img.w, height=img.h, dpi=img.dpi)
  }
  if (return.plot) { return(the.plot) }
}
