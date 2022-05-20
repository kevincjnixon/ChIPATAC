
addRank<-function(x, enrichment=FALSE){
  #Add a rank based on either the Log.P.Value (ascending)
  #Or calculate enrichment based on background and use that if enrichment=TRUE
  if(isTRUE(enrichment)){
    #Enrichment is %target with motif/%background with motif
    #Convert percentage columns to numeric
    t<-as.numeric(stringr::str_remove(x$X..of.Target.Sequences.with.Motif, "%"))
    b<-as.numeric(stringr::str_remove(x$X..of.Background.Sequences.with.Motif, "%"))
    x$enrichment<-t/b
    x<-x[order(x$enrichment, decreasing=T),]
    x$rank<-1:nrow(x)
  }else {
    x<-x[order(x$Log.P.value),]
    x$rank<-1:nrow(x)
  }
  return(x)
}


TF_plot<-function(x,y, title="test", xlab="x",ylab="y", lab=T, spec=NULL, numLab=5){
  #print(nrow(x))
  #Combine the dataframes based on ranks and calculate the differences:
  take<-c("Motif.Name","Log.P.value","rank")
  x<-x[,which(colnames(x) %in% take)]
  y<-y[,which(colnames(y) %in% take)]

  #merge
  res<-merge(x,y, by.x="Motif.Name",by.y="Motif.Name")
  #print(nrow(res))
  #Calculate the rank differences
  res$diff<-res$rank.y-res$rank.x
  res<-res[order(res$diff),]
  #Apply the colours
  cols<-colorRampPalette(c("blue","white","red"))(nrow(res))
  res$colour<-cols
  to_lab<-c()
  if(numLab>0){
    to_lab<-res[order(abs(res$diff), decreasing=T),][c(1:numLab),]
  }
  if(!is.null(spec)){
    for(i in 1:length(spec)){
      to_lab<-rbind(to_lab, res[grep(spec[i], res$Motif.Name, ignore.case=T),])
    }
  }
  #Now plot based on ranks
  par(mar=c(8,4,4,2), xpd=F)
  plot(res$rank.x, res$rank.y, pch=16, col=res$colour,  main=title, xlab=xlab, ylab=ylab)
  abline(a=0, b=1, lty=2)
  par(xpd=T)
  if(isTRUE(lab)){
    text(to_lab$rank.x, to_lab$rank.y, labels=sapply(strsplit(to_lab$Motif.Name, split="(", fixed=T),'[[',1), pos=4)
  }
  legend("bottom", inset=c(0,-0.5),
         legend=c(paste("Enriched in",xlab),paste("Enriched in",ylab)),
         pch=16, col=c("red","blue"))
  return(res[order(abs(res$diff), decreasing=T),])
}

#' Function to compare two sets of HOMER known motif enrichment results
#' Results from HOMER known TF motif enrichment analyses are ranked either by enrichment or significance and plotted for comparison
#' @param xfile string of path to knownResults.txt file of a HOMER motif enrichment analysis
#' @param yfile string of path to knownResults.txt file of another HOMER motif enrichment analysis
#' @param title String indicating title of figure
#' @param xlab String indicating the x-axis label (describing xfile). Default="x"
#' @param ylab String indicating the y-axis label (describing yfile). Default="y"
#' @param lab Boolean indicating if certain TFs should be labeled
#' @param spec String vector indicating TFs of interest that should be labeled
#' @param numLab Numeric indicating the number of top ranked TFs to be labeled (default=5)
#' @param enrichment Boolean indicating if the results should be ranked by enrichment. FALSE ranks results by significance. Default=FALSE.
#' @param returnRes Boolean indicating if data frame of results should be returned. Default=F.
#' @param return Plot comparing HOMER results. Data frame of TFs and corresponding rankings and differences between those rankings.
#' @export

compHOMER<-function(xfile, yfile, title="", xlab="x", ylab="y", lab=T, spec=NULL, numLab=5, enrichment=F, returnRes=F){
  message("Reading in files...")
  files<-c(xfile,yfile)
  TF_df<-list()
  for(i in 1:length(files)){
    if(file.exists(files[i])){
      TF_df[[i]]<-read.delim(files[i])
      #names(TF_df)[i]<-path1[i]
    }
  }
  to_remove<-c()
  for(i in 1:length(TF_df)){
    if(length(TF_df[[i]])==0){
      to_remove<-c(to_remove,i)
    }
  }
  if(length(to_remove)>0){
    TF_df<-TF_df[-to_remove]
  }
  if(length(TF_df)<2){
    stop("Only 1 file read in. Cannot plot without a second...")
  }
  if(length(TF_df)>2){
    message("More than 2 files read in. Using only first two.")
  }

  TF_rank<-lapply(TF_df, addRank, enrichment=enrichment)

  res<-TF_plot(TF_rank[[1]], TF_rank[[2]], title, xlab, ylab, lab, spec, numLab)
  if(isTRUE(returnRes)){
    return(res)
  }
}
