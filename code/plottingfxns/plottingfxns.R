 library(reshape2)
library(viridis)


colorvis <- function(COL){
  plot(NULL, xlim=c(0,length(COL)), ylim=c(0,1), 
       xlab="", ylab="", xaxt="n", yaxt="n")
  rect(0:(length(COL)-1), 0, 1:length(COL), 1, col=COL)
}

MatrixPrep <- function(v){
  # Melt a matrix and display upright, 1,1 in the corner
  numClusters = nrow(v)
  v <- melt(t(v)[,numClusters:1])
  v
}

MakeColormap <- function(v){
  make0 <- seq(min(v$value),max(v$value),length = 1000)
  make0 <- which(diff(sign(make0)) != 0)
  customcmap <- viridis(1000,option = 'plasma')
  customcmap[(make0-1):(make0+1)] <- "grey50"
  customcmap
}

fliplr <- function(x){
  x <- x[length(x):1]
  return(x)
}

p.signif.matrix <- function(p){
  # take matrix of p values and make asterisks
  #   ns: p > 0.05
  # *: p <= 0.05
  # **: p <= 0.01
  # ***: p <= 0.001
  # ****: p <= 0.000001
  p.new <- matrix(data = '',nrow = nrow(p),ncol=ncol(p),dimnames=dimnames(p))
  p.new[p > 0.05] <- ''
  p.new[p < 0.05 & p > 0.01] <- '*'
  p.new[p <= 0.01 & p > 0.001] <- '**'
  p.new[p <= 0.001 & p > 0.000001] <- '***'
  p.new[p <= 0.000001] <- '****'
  return(p.new)
  
}

imagesc <- function(X,caxis_name='',cmap='plasma',caxis_labels=NULL,clim=c(min(X,na.rm=T),max(X,na.rm=T)),
  xlabel='',ylabel='',yticklabels=rownames(X),xticklabels=as.character(colnames(X)),ttl='',noticks=FALSE,overlay = NULL,overlay.text.col='black',overlay.text.sz=2.5){
  # INPUTS:
  # X: matrix with dim names
  # cmap: name of colormap, for R colorbrewer
  # caxis_labels: tick labels for colorbar
  # clim: color axis limits
  # xlabel, ylabel, yticklabels, xticklabels, ttl: x axis labe, y axis label, y axis tick labels, x axis tick labels, plot title
  # noticks: logical option to remove all ticks
  # overlay: optional matrix of a text overlay (often p-value cut off labels or numbers), must be same size as X
  # OUTPUTS:
  # heatmap plot of X in style of matlab imagesc

  if(is.null(caxis_labels)){ # default axis label breaks
    caxis_breaks<-labeling::extended(clim[1], clim[2], m = 5)
    caxis_labels<-as.character(labeling::extended(clim[1], clim[2], m = 5))
  } else if(is.character(caxis_labels)){ # if axis is discrete then autogenerate breaks and label with provided labels
    caxis_breaks<-labeling::extended(clim[1], clim[2], m = length(caxis_labels))
  } else if(is.numeric(caxis_labels)){
    caxis_breaks<-labeling::extended(clim[1], clim[2], m = length(caxis_labels))
  }
  if(length(yticklabels) == 0){yticklabels <- rownames(X) <- as.character(1:nrow(X))}
  if(length(xticklabels) == 0){xticklabels <- colnames(X) <- as.character(1:ncol(X))}

  if(!is.null(overlay)){
    if(!identical(dim(overlay),dim(X))){print('ERROR: overlay values are different size than matrix'); return()}
  }
  

  X[X>max(clim)] <- max(clim) # threshold data based on color axis
  X[X < min(clim)] <- min(clim)

  melt_mat <- melt(t(X))
  melt_mat$Var2[is.na(melt_mat$Var2)] <- 'NA'
  melt_mat$Var1 <- as.character(melt_mat$Var1)
  p<-ggplot() + geom_tile(data = melt_mat, aes(x=Var1, y=Var2, fill=value)) + 
    scale_x_discrete(limits=as.character(colnames(X)),labels=xticklabels,expand = c(0,0)) +
    scale_y_discrete(limits=rev(rownames(X)),labels=rev(yticklabels),expand = c(0,0))
  if(cmap =='plasma'){
    p <- p + scale_fill_viridis(option = 'plasma',name=caxis_name,limits=clim,breaks=caxis_breaks,labels=caxis_labels)
  } else if(cmap == 'redblue'){
    p <- p + scale_fill_gradientn(colours = c('#8B0000','#c23b22','#ffffff','#779ecb','#00008b'),
                           guide = "colorbar", limits=clim,breaks=caxis_breaks,
                           na.value = 'grey',name=caxis_name)
  } else if(cmap == 'redblue_asymmetric'){
    p <- p + scale_fill_gradientn(colours = c('#8B0000','#c23b22','#ffffff','#779ecb','#00008b'),
                           values = scales::rescale(c(clim[1],0,clim[2])),
                           guide = "colorbar", limits=clim,breaks=caxis_breaks,
                           na.value = 'white',name=caxis_name)
  } else if(cmap == 'custom1_asymmetric'){
    p <- p + scale_fill_gradientn(colours = c('#c87157','#ffffff','#57c8bd'),
                                  values = scales::rescale(c(clim[1],0,clim[2])),
                                  guide = "colorbar", limits=clim,breaks=caxis_breaks,
                                  na.value = 'white',name=caxis_name)
  } else {
    pal.idx <- which(rownames(brewer.pal.info) == cmap)  
    cols <- brewer.pal(brewer.pal.info$maxcolors[pal.idx], cmap)
    p <- p + scale_fill_gradientn(colours = cols,
                           guide = "colorbar", limits=clim,
                           na.value = 'grey',name=caxis_name,breaks=caxis_breaks,labels=caxis_labels)
  } 
  p <- p + theme_bw()
  p <- p + xlab(xlabel)+ylab(ylabel)+ ggtitle(ttl)+
      theme(text=element_text(size=8),plot.title=element_text(hjust=0.5,size=8,face='bold'))
  if(noticks){p <- p + theme_void()}
  if(!is.null(overlay)){
    melt_ov_mat <- melt(t(overlay))
    melt_ov_mat$Var1 <- as.character(melt_ov_mat$Var1)
    p <- p + geom_text(data = melt_ov_mat, aes(x=Var1, y=Var2, label=value),size=overlay.text.sz,color=overlay.text.col)
  }
  return(p)

}

nice_cbar <- function(pos='right'){
  if(pos=='bottom' | pos == 'top'){
    thm <- theme(legend.box = 'none',legend.background = element_blank(),
            legend.margin = ggplot2::margin(0,0,0,0),
            legend.position = pos,legend.key.height = unit(0.1,units='cm'),legend.key.width=unit(1.25,'cm'))
  } else if(pos=='left' | pos == 'right'){
    thm <- theme(legend.box = 'none',legend.background = element_blank(),
            legend.margin = ggplot2::margin(0,0,0,0),
            legend.position = pos,legend.key.width=unit(0.1,'cm'))
  }
  return(thm)
}
p.violin.nptest <- function(val,dist,xlab='',ylab='',col='blue'){
  # val: observed value
  # dist: null distribution
  # xlab: name of quantity for x-axis tick label

  # plots violin plot of distribution, with val as a point
  df <- data.frame(x=rep(xlab,length(dist)),y=dist)
  p <- ggplot(df) + geom_violin(aes(x=x,y=dist),fill=col,alpha=0.5) + geom_point(aes(x=xlab,y=val)) +
  xlab('') + theme_classic() + theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) + theme(text=element_text(size=8))
  return(p)
}


getGroupColors <- function(grps){
  # return colors to use for PNC and 22q respectively
  # return(c('#005C9F','#FF8400')) # blue orange from n-back paper
   return(c('#62589F','#93AD90'))
  #return(c('#63ae61','#683388'))
  #return(c('#b89a3e','#683388'))
}

p.group.jitter <- function(df,xname,yname,grpname,
                           ylabel='',xlabel='',
                           xticklabels=unique(df[,xname]),cols=rep('k',length(unique(df[,grpname]))),
                           alpha=1,ylim=range(df[,yname],na.rm = T)){
  # INPUTS:
  # df: dataframe with 3 columns for measured values, the type of value, and groups
  # xname: character column name in df containing type of values
  # yname: character column name in df containing values
  # grpname: character column name in df containing groups
  # i.e. FO (yname) for 5 clusters (xname) for PNC and 22q (grpname)
  # other inputs are ancillary plotting visuals
  
  # OUTPUTS:
  # p: ggplot object
  p <- ggplot(data=df) + geom_jitter(aes(x = eval(parse(text=xname)),y = eval(parse(text=yname)),color = eval(parse(text=grpname))),stroke=0,alpha=alpha,position = position_jitterdodge()) + theme_classic() +
    scale_color_manual(limits = rev(unique(df[,grpname])), values = cols) + 
    ylab(ylabel) + xlab(xlabel) +theme(text = element_text(size = 8)) +
    theme(legend.title = element_blank()) + 
    scale_y_continuous(limits = ylim)+
    scale_x_discrete(limits = as.character(unique(df[,xname])), breaks=as.character(unique(df[,xname])), labels = xticklabels) +
    theme(legend.key.size = unit(0.5,'line'))
  return(p)
}

p.single.jitter <- function(df,yname,grpname,
                           ylabel='',xlabel='',xticklabels=unique(df[,grpname]),
                           cols=rep('k',length(unique(df[,grpname]))),
                           alpha=1,ylim=range(df[,yname],na.rm = T)){
  # INPUTS:
  # df: dataframe with 3 columns for measured values, the type of value, and groups
  # xname: character column name in df containing type of values
  # yname: character column name in df containing values
  # grpname: character column name in df containing groups
  # i.e. FO (yname) for 5 clusters (xname) for PNC and 22q (grpname)
  # other inputs are ancillary plotting visuals
  
  # OUTPUTS:
  # p: ggplot object
  p <- ggplot(data=df) + geom_jitter(aes(x = eval(parse(text=grpname)),y = eval(parse(text=yname)),color = eval(parse(text=grpname))),stroke=0,alpha=alpha,position = position_jitterdodge()) + theme_classic() +
    scale_color_manual(limits = rev(unique(df[,grpname])), values = cols) + 
    ylab(ylabel) + xlab(xlabel) +theme(text = element_text(size = 8)) +
    theme(legend.title = element_blank()) + 
    scale_y_continuous(limits = ylim)+
    scale_x_discrete(limits = as.character(unique(df[,grpname])), breaks=as.character(unique(df[,grpname])), labels = xticklabels) +
    theme(legend.key.size = unit(0.5,'line'))
  return(p)
}

p.xy <- function(x,y,xlab,ylab,ttl='',col='black'){
  df <- data.frame(x=x,y=y)
  r.text <- paste('r_p =',signif(cor(x,y,use='pairwise.complete.obs'),2),'r_s =',signif(cor(x,y,method='spearman',use='pairwise.complete.obs'),2))
  p <- ggplot(df) + geom_point(aes(x=x,y=y),color=col,alpha=0.3,size=0.8,stroke=0) + geom_smooth(aes(x=x,y=y),fill=col,color=col,method='lm') + 
    xlab(xlab) + ylab(ylab) + ggtitle(ttl) + 
    annotate("text",size = 2, x = min(df$x),y = min(df$y,na.rm=TRUE), label = r.text,hjust=0) +
      theme_classic() + theme(text = element_text(size = 8)) + 
      theme(plot.title = element_text(size=8,hjust=0.5,face = "bold")) +
      theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) + theme(legend.position = 'none')
  return(p)
}

p.xy.flex <- function(x,y,xlab,ylab,ttl='',col='black',alpha=1,r.method='pearson',sm.method = 'lm',formula=NULL,
                      p=NULL,r=NULL,ptxt='p = ',rtxt='r = ',parse=F,pt.sz=1,ln.sz=1,pos = 'bottom.right'){
  # INPUTS:
  # x: x variable, vector
  # y: y variable, vector
  # xlab, ylab, ttl: character labels for axes
  # col: point color and line color, RGB hex code or R native color or RGB value
  # alpha: opacity of points
  # sm.method: method for geom_smooth to draw line of best fit
  # formula: if using gam or nonlinear fit, provide formula\
  # p,r: p-value or stat/coefficient to display, for instance, if you post-hoc correct outside the function
  # ptxt,rtxt: what to say for p=, i.e. p[FDR]== or just p=
  # parse: logical specifying whether to parse text for annotations. if parsing should use 'p ==' instead of 'p = '
  # pos: (top/bottom).(right/left) character i.e. 'top.right' or 'bottom.left' specifying location of label
  #
  # OUTPUTS:
  # scatter plot with r and p value for pearson correlation between x and y
  # and linear fit

  df <- data.frame(x=x,y=y)
  c.test <- cor.test(df$x,df$y,method=r.method,use='pairwise.complete.obs')
  if(is.null(p)){p <- c.test$p.value}
  if(is.null(r)){r <- c.test$estimate}
  r.text <- paste0(rtxt,signif(r,2))
  p.text <- paste0(ptxt,signif(p,2)) # annotation
  pos <- strsplit(pos,'[.]')[[1]] # get vertical (v) and horizontal (h) position
  v <- pos[1]; h <- pos[2] 
  if(v == 'top'){v <- 1;scl.r<-1;scl.p<-2}else{v <- -1;scl.r<-2;scl.p<-1}
  if(h == 'right'){h <- 1}else{h <- -1}

  p <- ggplot(df) + geom_point(aes(x=x,y=y),color=col,stroke=0,alpha=alpha,size=pt.sz) + 
  geom_smooth(aes(x=x,y=y),fill=col,color=col,method=sm.method,formula=formula,size=ln.sz) + 
    xlab(xlab) + ylab(ylab) + ggtitle(ttl) + 
    annotate("text",size = 2, x = h*Inf,y =v*Inf, label = r.text,hjust=1*h,vjust=scl.r*v,parse=parse) +
    annotate("text",size = 2, x = h*Inf,y =v*Inf, label = p.text,hjust=1*h,vjust=scl.p*v,parse=parse) +
      theme_classic() + theme(text = element_text(size = 8)) + 
      theme(plot.title = element_text(size=8,hjust=0.5,face = "bold")) +
      theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) + theme(legend.position = 'none')
  return(p)
}

