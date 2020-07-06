data.summary=function(x,plot=TRUE)
{
  require(car)
  x=as.vector(x)
  na.count=sum(is.na(x))
  if(na.count>0){
    x=na.omit(x)
    if(plot==TRUE)
    {
      cat("A stem and leaf diagram is","\n")
      stem(x)
    }
    cat("The number of missing data is:", na.count,"\n")
    cat("The summary below is obtained after deleting the missing data","\n")
    cat("\n")
    }else{
      if(plot==TRUE)
      {
        cat("A stem and leaf plot is","\n")
        stem(x)
      }
      cat("There is no missing value.","\n")
      cat("\n")
    }
  cat("Summary:","\n")
  res=c(min(x),mean(x),var(x),sd(x),max(x),max(x)-min(x),as.numeric(quantile(x,0.25)),median(x),as.numeric(quantile(x,0.75)),as.numeric(quantile(x,0.75)-quantile(x,0.25)))
  result=round(res,4)
  statistics=c("min","mean","variance","std","max","range","Q1","Median","Q3","IRQ")
  print.table(rbind(statistics,result))
  cat("\n")
  if(plot==TRUE)
  {
    par(mfrow=c(2,2))
    par(mar=c(4,4,2,.1))
    plot(x,main="Scatterplot of x")
    boxplot(x,main="boxplot of x",ylab="x")
    hist(x)
    qqPlot(x,xlab="normal quantile",ylab="Observations")
    par(mfrow=c(1,1))
  }
}
