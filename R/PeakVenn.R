getOL_gr<-function(x, retVals=F){
  if(length(x)>5){
    stop("max length of x is 5...")
  }
  OL=list()
  index<-1
  #start with 1 on one comparisons:
  for(i in 1:length(x)){
    if(i < length(x)){
      for(j in (i+1):length(x)){
        #OL[[index]]<-x[[i]][which(x[[i]] %in% x[[j]])]
        OL[[index]]<-subsetByOverlaps(x[[i]],x[[j]])
        names(OL)[index]<-paste0("n",i,j)
        index<-index+1
      }
    }
  }
  if(length(x)>2){
    #Now three way comparisons
    for(i in 1:length(x)){
      if(i<length(x)){
        for(j in (i+1):length(x)){
          if(j<length(x)){
            #tmp<-x[[i]][which(x[[i]] %in% x[[j]])]
            tmp<-subsetByOverlaps(x[[i]],x[[j]])
            for(k in (j+1):length(x)){
              #OL[[index]]<-tmp[which(tmp %in% x[[k]])]
              OL[[index]]<-subsetByOverlaps(tmp,x[[k]])
              names(OL)[index]<-paste0("n",i,j,k)
              index<-index+1
            }
          }
        }
      }
    }
  }
  if(length(x)>3){
    #Now for the four way comparisons
    for(i in 1:length(x)){
      if(i<length(x)){
        for(j in (i+1):length(x)){
          if(j<length(x)){
            #tmp<-x[[i]][which(x[[i]] %in% x[[j]])]
            tmp<-subsetByOverlaps(x[[i]],x[[j]])
            for(k in (j+1):length(x)){
              if(k<length(x)){
                #tmp<-tmp[which(tmp %in% x[[k]])]
                tmp<-subsetByOverlaps(tmp, x[[k]])
                for(d in (k+1):length(x)){
                  #OL[[index]]<-tmp[which(tmp %in% x[[d]])]
                  OL[[index]]<-subsetByOverlaps(tmp,x[[d]])
                  names(OL)[index]<-paste0("n",i,j,k,d)
                  index<-index+1
                }
              }
            }
          }
        }
      }
    }
  }
  if(length(x)==5){
    #Now for the final comparison:
    tmp<-x[[1]]
    for(i in 2:length(x)){
      #tmp<-tmp[which(tmp %in% x[[i]])]
      tmp<-subsetByOverlaps(tmp, x[[i]])
    }
    OL[[index]]<-tmp
    names(OL)[index]<-"n12345"
  }
  #print(names(OL))
  if(isFALSE(retVals)){
    OL<-as.list(lengths(OL))
  }
  return(OL)
}

makeNames<-function(x){
  y<-c()
  for(i in 1:length(x)){
    if(i < length(x)){
      for(j in (i+1):length(x)){
        y<-c(y, paste(x[i],x[j],sep="_"))
      }
    }
  }
  if(length(x)>2){
    for(i in 1:length(x)){
      if(i<length(x)){
        for(j in (i+1):length(x)){
          if(j<length(x)){
            for(k in (j+1):length(x)){
              y<-c(y, paste(x[i],x[j],x[k], sep="_"))
            }
          }
        }
      }
    }
  }
  if(length(x)>3){
    for(i in 1:length(x)){
      if(i<length(x)){
        for(j in (i+1):length(x)){
          if(j<length(x)){
            for(k in (j+1):length(x)){
              if(k<length(x)){
                for(d in (k+1):length(x)){
                  y<-c(y, paste(x[i],x[j],x[k],x[d],sep="_"))
                }
              }
            }
          }
        }
      }
    }
  }
  if(length(x)==5){
    y<-c(y, paste(x[1],x[2],x[3],x[4],x[5], sep="_"))
  }
  return(y)
}

#' Function to plot Venn diagram of overlapping peaks
#' @param x named list of GenomicRagnes objects to compare. Maximum 5 sets of regions.
#' @param title String indicating title of plot
#' @param cols colour palette to use. Can be RColorBrewer Palette, or string of rgb, hexadecimal, or colour names the same length of 'x'
#' @param lty Line type. Default="blank" for no lines. 1 for solid lines, 2 for dashed lines, etc.
#' @param scale Boolean indicating if circles should be scaled to the size of each region. Works for comaparisons of less than 4.
#' @param retVals Boolean indicating if overlap regions should be returned. Default=FALSE.
#' @return Venn diagram plot and if retVals=T, a named list of GenomicRanges for each overlapping region
#' @export

plotVenn_gr<-function(x, title="Venn Diagram", cols="Dark2", lty="blank", scale=F, retVals=F){
  require(VennDiagram, quietly=T)
  y<-getOL_gr(x)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(2,1, heights=unit(c(0.25, 10),"null"))))
  venn<-NULL
  if(length(x)==1){
    venn<-draw.single.venn(area=length(x),
                           fill=BinfTools::colPal(cols)[1:length(x)],
                           lty=lty,
                           category=names(x))
  }
  if(length(x)==2){
    venn<-draw.pairwise.venn(area1=length(x[[1]]),
                             area2=length(x[[2]]),
                             cross.area=y$n12,
                             fill=BinfTools::colPal(cols)[1:length(x)],
                             lty=lty,
                             category = names(x),
                             scaled=scale,
                             euler.d=scale)
  }
  if(length(x)==3){
    if(isTRUE(scale)){
      assign("overrideTriple", TRUE, envir=.GlobalEnv)
    }
    venn<-draw.triple.venn(area1=length(x[[1]]),
                           area2=length(x[[2]]),
                           area3=length(x[[3]]),
                           n12=y$n12,
                           n23=y$n23,
                           n13=y$n13,
                           n123=y$n123,
                           fill=BinfTools::colPal(cols)[1:length(x)],
                           lty=lty,
                           category = names(x),
                           scaled=scale,
                           euler.d=scale)
    if(isTRUE(scale)){
      rm(overrideTriple, pos=".GlobalEnv")
    }
  }
  if(length(x)==4){
    venn<-draw.quad.venn(area1=length(x[[1]]),
                         area2=length(x[[2]]),
                         area3=length(x[[3]]),
                         area4=length(x[[4]]),
                         n12=y$n12, n13=y$n13, n14=y$n14,
                         n23=y$n23, n24=y$n24, n34=y$n34,
                         n123=y$n123, n124=y$n124, n134=y$n134, n234=y$n234,
                         n1234=y$n1234,
                         fill=BinfTools::colPal(cols)[1:length(x)],
                         lty=lty,
                         category=names(x))
  }
  if(length(x)==5){
    venn<-draw.quintuple.venn(area1=length(x[[1]]), area2=length(x[[2]]),
                              area3=length(x[[3]]), area4=length(x[[4]]), area5=length(x[[5]]),
                              n12=y$n12, n13=y$n13, n14=y$n14, n15=y$n15, n23=y$n23, n24=y$n24,
                              n25=y$n25, n34=y$n34, n35=y$n35, n45=y$n45, n123=y$n123, n124=y$n124,
                              n125=y$n125, n134=y$n134, n135=y$n135, n145=y$n145, n234=y$n234, n235=y$n235,
                              n245=y$n245, n345=y$n345, n1234=y$n1234, n1235=y$n1235, n1245=y$n1245,
                              n1345=y$n1345, n2345=y$n2345, n12345=y$n12345,
                              fill=BinfTools::colPal(cols)[1:length(x)],
                              lty=lty,
                              category=names(x))
  }
  #print(venn, vp=viewport(layout.pos.row=2))
  grid.text(title, vp=viewport(layout.pos.row=1))

  if(isTRUE(retVals)){
    res<-getOL_gr(x, retVals=retVals)
    #res<-c(nonOL(x, res), res)
    names(res)<-makeNames(names(x))
    return(res)
  }
}
