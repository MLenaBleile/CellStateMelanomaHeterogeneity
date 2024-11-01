library(ggplot2)
library(Seurat)
ScSpatialFeaturePlot = function(spat, features,slot="counts",cells=NULL,ncol=NULL,
                                facet_labeller=NULL,
                                logs=F,pt.size.factor=1.6, LegLabel="Count", uq=.999, lq=0, flip=1){
  fnames=features
  MinMaxq = function(X,uq=.99, lq=0.05){MinMax(X, min=quantile(X,lq), max=quantile(X,uq))}
  rawcoords=GetTissueCoordinates(spat)
  if(is.null(cells)){
    cells=Cells(spat)
  }
  if(flip==1){
    rawcoords$x = -rawcoords$imagecol
    rawcoords$y = rawcoords$imagerow
  }else if(flip==2){
    rawcoords$x = rawcoords$imagerow
    rawcoords$y = rawcoords$imagecol
  }else{
    rawcoords$x = rawcoords$imagecol
    rawcoords$y = -rawcoords$imagerow
  }
  all.plotdf=data.frame()
  for(fname in fnames){
    if(fname %in% colnames(spat@meta.data)){
      spat$intfeat = MinMaxq(spat[[fname]][,1], uq=uq, lq=lq)
                             }else{
        feat = FetchData(spat, slot=slot, vars=fname)[,1]
        spat$intfeat=MinMaxq(feat, uq=uq, lq=lq)
      }
    plotdf= data.frame(x=rawcoords$x, 
                       y= rawcoords$y,
                       feature=spat$intfeat[rownames(rawcoords)],
                       fname = rep(fname, nrow(rawcoords)),
                       #ftype = rep(names(fname), nrow(rawcoords)),
                       #fcols = color.gradient(spat$intfeat[rownames(rawcoords)],colors=red_grad, colsteps=10000),
                       cells=factor(Cells(spat)))
    all.plotdf = rbind(all.plotdf, plotdf)
  }

  if(logs){
    all.plotdf$feature=log(all.plotdf$feature+1)

  }
  
  #rownames(all.plotdf) = all.plotdf$cells
  all.plotdf = all.plotdf[all.plotdf$cells %in% cells,]
  p1=ggplot(all.plotdf)+geom_point(mapping=aes(x=x, y=y, color=feature), size=pt.size.factor)+
    #scale_color_manual(values=plotdf$scfred)+
    scale_color_viridis_c(option="turbo")+
    theme_bw()+
    
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+labs(x="", y="",color=LegLabel )+facet_wrap(~fname, ncol=ncol)
  return(p1)
}