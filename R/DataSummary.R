x=read.csv("https://raw.githubusercontent.com/Harrindy/StatEngine/master/Data/CompressiveStrength.csv")
data.summary=function(x)
{
  x=as.vector(x)
  na.count=sum(is.na(x))
  if(na.count>0){
    cat("The number of missing data is:", na.count)
    cat("The summary below is obtained after deleting the missing data")
    x=na.omit(x)
    }else{
      cat("There is no missing value.")
    }
  header=c("min","mean","variance","std","max","range")
  res=c(min(x),mean(x),var(x),sd(x),max(x),max(x)-min(x))
  cat(rbind(header,res))
  robu.header=c("Q1","Median","Q3","IRQ")
}

data.summary(x)
