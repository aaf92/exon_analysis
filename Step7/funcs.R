library(DESeq2)
library(dplyr)
library(tibble)

tscale=function(x){
  s = apply(x, 1, sd)
  m = apply(x, 1, mean)
  x = sweep(x, 1, m)
  x = sweep(x, 1, s, "/")
  x
}

tyscale=function(x){
  s = apply(x, 2, sd)
  m = apply(x, 2, mean)
  x = sweep(x, 2, m)
  x = sweep(x, 2, s, "/")
  x
}

splitNamesByToken=function(cnames,char="_"){
  nn=length((strsplit(cnames[1], char))[[1]])
  show(nn)
  matrix(unlist(strsplit(cnames, char)), ncol=nn, byrow=T)
}

getSizeFactors=function(counts){
  colnames(counts)=NULL
  tmp=DESeqDataSetFromMatrix(counts, design=~1,colData = data.frame(rep(1, ncol(counts))))
  tmp=estimateSizeFactors(tmp)
  sf=sizeFactors(tmp)
}

plotPCAggplot=function(data, grp=NULL, grp2=NULL,scale=T, pc1=1, pc2=2, mycol=NULL, rev1=F, rev2=F, sz=NULL, bottomright=NULL, return.pca=F, subset=NULL, text=NULL, svdres=NULL, display.subset=NULL, alpha=1){
  #  message(paste0("scale is set to ",scale))
  if(is.null(grp)){
    grp=rep("1", ncol(data))
  }
  if(is.null(colnames(data))){
    colnames(data)=1:ncol(data)
  }
  if(is.null(sz)){
    sz=3
    if(ncol(data)>200){
      sz=2
    }
    if(ncol(data)>1000){
      sz=1
    }
  }
  if(!is.null(subset)){
    data=data[,subset]
    
    
    grp=grp[subset]
    grp2=grp2[subset]
    if(!is.null(text)){
      text=text[subset]
    }
  }
  if(scale && is.null(svdres)){
    vv=apply(data,1, var)
    data=tscale(data[vv>1e-6,])
  }
  if (is.null(svdres)){
    if(ncol(data)>50){
      set.seed(123); svdres=rsvd(data, k=max(c(10,pc1,pc2)))
    }
    else{
      svdres=svd(data)
    }
  }
  tot=sum(svdres$d^2)
  # show(c(pc1, pc2))
  #  colnames(svdres$v)=paste0("PC", 1:ncol(svdres$v))
  if(!is.null(bottomright)){
    ii=which(grp==bottomright)
    if(length(ii)==0 & !is.null(grp2)){
      ii=which(grp2==bottomright)
    }
    if(mean(svdres$v[ii,pc1])<0)
      rev1=T
    if(mean(svdres$v[ii,pc2])>0)
      rev2=T
  }
  if(rev1){
    svdres$v[,pc1]=-svdres$v[,pc1]
  }
  if(rev2){
    svdres$v[,pc2]=-svdres$v[,pc2] 
  }
  if(!return.pca){
    vare=(svdres$d^2)/tot
    #show(svdres$d)
    if(is.null(mycol)){
      mycol=c(brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(12,"Set3"))
      mycol[6]="black"
    }
    if(is.null(grp2)){
      tmpall=tmp=list(PC1=svdres$v[,pc1], PC2=svdres$v[,pc2], grp=(grp), names=colnames(data))
      tmp=as.data.frame(tmp)
      if(!is.null(text)){
        tmp$text=text
      }
      if(!is.null(display.subset)){
        tmp=tmp[display.subset,]
        tmp[-display.subset,"PC1"]=NA
      }
      p=ggplot()+geom_point(data=tmp, aes(x=PC1, y=PC2, color=grp), size=sz, alpha=alpha)+theme_bw()+labs(x=paste0("PC", pc1, "(", round(100*vare[pc1],2), "%)"), y=paste0("PC",pc2, "(", round(100*vare[pc2],2), "%)"))
      # show(scale_type(tmp$grp))
      
      if(scale_type(tmp$grp)=="discrete"){
        p=p+scale_colour_manual(values=mycol)
      }
    }
    else{
      tmpall=tmp=list(PC1=svdres$v[,pc1], PC2=svdres$v[,pc2], grp=(grp),grp2=(grp2),names=colnames(data))
      tmp=as.data.frame(tmp)
      if(!is.null(text)){
        tmp$text=text
      }
      if(!is.null(display.subset)){
        #tmp=tmp[display.subset,]
        
        tmp[-display.subset,"PC1"]=NA
      }
      p=ggplot()+geom_point(data=tmp, aes(x=PC1, y=PC2, color=grp, shape=grp2), size=sz, alpha=alpha)+theme_bw()+scale_colour_manual(values=mycol)+scale_shape_manual(values=c(19,15, 17,10,9, 8,7,6,1,2,3))+labs(x=paste0("PC", pc1, "(", round(100*vare[pc1],2), "%)"), y=paste0("PC", pc2, "(", round(100*vare[pc2],2), "%)"))
      #+
      #scale_shape(solid=F)
    }
    if(!is.null(text)){
      p+geom_text_repel(data=tmp, aes(x=PC1, y=PC2, label=text), size=10, box.padding = 1)
    }
    else{
      p+xlim(range(tmpall$PC1))+ylim(range(tmpall$PC2))
    }
  }
  else{
    invisible(svdres)
  }
}

normalizeToControls=function(data, batch, isControl){
  myGSEs=list()
  for (g in unique(batch)){
    ii=which(batch==g)
    datasub=data[, ii]
    myGSEs[[g]]$genes=data[,ii]
    myGSEs[[g]]$pheno=as.data.frame(list(isNotControl=as.numeric(!isControl[ii])))
    rownames(myGSEs[[g]]$pheno)=colnames(data)[ii]
    myGSEs[[g]]$pheno$platform="seq"
  }
  
  cres=COCONUT(GSEs = myGSEs, control.0.col = "isNotControl", byPlatform = FALSE)
  #the output splits controls and other so we put it back together
  dataC <- as.matrix( Reduce(cbind, lapply(cres$COCONUTList, function(x) x$genes)))
  dataC2 <- as.matrix( Reduce(cbind, lapply(cres$controlList$GSEs, function(x) x$genes)))
  #combine
  dataC=cbind(dataC, dataC2)
  #put back in the original order
  dataC=dataC[, colnames(data)]
  
}

MAplot <- function(dir, compare, res) {
  fname <- paste0(dir, "/MA_plot_", compare, ".pdf")
  pdf(fname)
  plotMA(res, ylim=c(-3,3))
  dev.off()
}

writeResults <- function(dir, compare, res) {
  write.csv(as.data.frame(resOrdered), 
            file=paste0(dir, "/resultsOrdered_", compare, ".txt"))
}

svd_filter <- function(y, x, ...) {
  # pca to get top 10 percent features from pc1
  svdres=svd(tyscale(x))
  V <- svdres$v
  first_pc_loadings <- V[, 1]
  names(first_pc_loadings) <- colnames(x)

  # Sort the gene loadings in decreasing order
  sorted_loadings <- sort(first_pc_loadings, decreasing = TRUE, key = function(x) abs(x))
  order_indices <- order(abs(first_pc_loadings), decreasing = TRUE)
  sorted_loadings <- first_pc_loadings[order_indices]


  num_elements <- length(sorted_loadings)
  first_10_percent <- sorted_loadings[1:ceiling(0.1 * num_elements)]

  first_10_percentNames <- names(first_10_percent)
  first_10_percentNames

  # return column indices for features to retain
  column_indices <- which(names(x) %in% first_10_percentNames)
  column_indices

}

filter_x <- function(y, x, percent=0.1, ...) {
  svdres=svd(tyscale(x))
  pc1 <- svdres$v[,1]
  abs_pc1 <- abs(pc1)
  sortedPC1_indices <- order(abs_pc1, decreasing = TRUE)
  top_10_percent_index <- round(length(sortedPC1_indices) * percent)
  top_10_percent_indices <- sortedPC1_indices[1:top_10_percent_index]
  top_10_percent_indices
}