#' Differential bursting analysis
#'
#' This function perfoms the differential bursting analysis of bursting frequency(i.e. kon) and bursting size(s/koff)
#' between two study conditions (e.g. Tumor vs. Normal) using
#' MCMC output from BISC_estimate().
#'
#' @export
#' @param data1 output from BISC_estimate() of treatment condition.
#' @param data2 output from BISC_estimate() of control condition.
#' @param frequency testing the H0: kon=frequency.
#' @param size testing the H0: size=size.
#' @param log2 compute difference of bursting estimates in log2 scale: T or F
#' @return This function outputs a list object contains testing results for three parameters (i.e. kon, koff and s).
#'
#'

DB=function(data1,data2,frequency,size,log2){
  test1=data1$MCMC
  test2=data2$MCMC
  N=dim(data$estimation)[1]
  if(log2==T){
    kondif=log2(test1[,1:N]/test2[,1:N])
    sizedif=log2(test1[,(2*N+1):3*N]/test1[,(N+1):2*N]/(test2[,(2*N+1):3*N])/test2[,(N+1):2*N])
  }else{
    kondif=test1[,1:N]-test2[,1:N]
    sizedif=test1[,(2*N+1):3*N]/test1[,(N+1):2*N]-(test2[,(2*N+1):3*N])/test2[,(N+1):2*N])

  }

  test_frequency=matrix(data = NA,nrow = N,ncol=3)
    for (i in 1:N) {
      test_frequency[i,]=c(mean(kondif[,i]<kona),mean(kondif[,i]>kona),mean(kondif[,i]))
    }
    colnames(test_kon)=c("Ha:kon>kona","Ha:kon<kona","log2FC")

    test_size=matrix(data = NA,nrow = N,ncol=3)
    for (i in 1:N) {
      test_size[i,]=c(mean(sizedif[,i]<kona),mean(sizedif[,i]>kona),mean(sizedif[,i]))
    }
    colnames(test_size)=c("Ha:size>sizea","Ha:size<sizea","log2FC")

  return(list(test_frequency=test_frequency,test_size=test_size))
}

volcano_plot=function(file,comparison,name){
  pdf(file = file)
  volcanoData <- data.frame(logFC=as.numeric(comparison$table$logFC), negLogPval=-log10(comparison$table$FDR))

  plot(volcanoData, pch=20,col="grey",ylab="q valuev(-log10)",
       xlab = "Fold change (log2)",main=name)
  abline(h=1.30103,lty=2)
  axis(labels = "q<0.05",side = 2, at = 1.30103, las = 2,padj = 0)
  pos=volcanoData[which(volcanoData$logFC>1 & volcanoData$negLogPval>1.30103),]
  neg=volcanoData[which(volcanoData$logFC< -1 & volcanoData$negLogPval>1.30103),]
  points(neg, pch=20, col="blue")
  points(pos, pch=20, col="red")
  legend("topleft", pch=c(20,20), col=c("red", "blue"), c(paste(dim(pos)[1]," upregulated genes",sep=""),paste(dim(neg)[1], " downregulated genes",sep="")),  box.col="darkgreen", cex=1.2)

}
