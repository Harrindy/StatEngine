Zinterval.data=function(x,sigma,level=0.95)
{
  if(level>=1|level<0){return("the confidence level should be between 0 and 1!")}
  alpha=1-level
  z_a=qnorm(1-alpha/2)
  n=length(x)
  barx=mean(x)
  cat("The sample mean is", mean(x),"and sample size is", n,"\n")
  cat("A", level*100, "% two-sided confidence interval for the population mean is [", barx-z_a*sigma/sqrt(n),",",barx+z_a*sigma/sqrt(n),"]","\n")
  z_a=qnorm(1-alpha)
  cat("A", level*100, "% upper-confidence bound for the population mean is", barx+z_a*sigma/sqrt(n),"\n")
  cat("A", level*100, "% lower-confidence bound for the population mean is", barx-z_a*sigma/sqrt(n),"\n")
}

Zinterval.stat=function(n,barx,sigma,level=0.95)
{
  if(level>=1|level<0){return("the confidence level should be between 0 and 1!")}
  alpha=1-level
  z_a=qnorm(1-alpha/2)
  cat("The sample mean is", barx,"and sample size is", n,"\n")
  cat("A", level*100, "% two-sided confidence interval for the population mean is [", barx-z_a*sigma/sqrt(n),",",barx+z_a*sigma/sqrt(n),"]","\n")
  z_a=qnorm(1-alpha)
  cat("A", level*100, "% upper-confidence bound for the population mean is", barx+z_a*sigma/sqrt(n),"\n")
  cat("A", level*100, "% lower-confidence bound for the population mean is", barx-z_a*sigma/sqrt(n),"\n")
}

sample.size.Zinterval=function(E,sigma,level=0.95)
{
  if(level>=1|level<0){return("the confidence level should be between 0 and 1!")}
  alpha=1-level
  z_a=qnorm(1-alpha/2)
  n=ceiling((z_a*sigma/E)^2)
  cat("The samllest sample size to control the margin of error E <=", E, "at",level*100, "% confidence level is", n, "\n")
}

AZinterval.data=function(x,level=0.95)
{
  s=sd(x)
  if(level>=1|level<0){return("the confidence level should be between 0 and 1!")}
  alpha=1-level
  z_a=qnorm(1-alpha/2)
  n=length(x)
  barx=mean(x)
  cat("The sample mean is", mean(x),", sample standard deviation is", s, ", and sample size is", n,"\n")
  cat("A large-sample confidence interval for the population mean with confidence level of approximately", level*100, "% is [", barx-z_a*s/sqrt(n),",",mean(x)+z_a*sigma/sqrt(n),"]","\n")
  z_a=qnorm(1-alpha)
  cat("A large-sample upper-confidence bound for the population mean with confidence level of approximately", level*100, "% is", barx+z_a*s/sqrt(n),"\n")
  cat("A large-sample lower-confidence bound for the population mean with confidence level of approximately", level*100, "% is", barx-z_a*s/sqrt(n),"\n")
}

AZinterval.stat=function(n,barx,s,level=0.95)
{
  if(level>=1|level<0){return("the confidence level should be between 0 and 1!")}
  alpha=1-level
  z_a=qnorm(1-alpha/2)
  cat("The sample mean is", mean(x),", sample standard deviation is", s, ", and sample size is", n,"\n")
  cat("A large-sample confidence interval for the population mean with confidence level of approximately", level*100, "% is [", barx-z_a*s/sqrt(n),",",mean(x)+z_a*sigma/sqrt(n),"]","\n")
  z_a=qnorm(1-alpha)
  cat("A large-sample upper-confidence bound for the population mean with confidence level of approximately", level*100, "% is", barx+z_a*s/sqrt(n),"\n")
  cat("A large-sample lower-confidence bound for the population mean with confidence level of approximately", level*100, "% is", barx-z_a*s/sqrt(n),"\n")
}

t.quantile=function(df,q)
{
  if(df<=0){return("df must be positive")}
  if(q<=0|q>=1){return("q must be between 0 and 1")}
  return(qt(q,df))
}


Tinterval.data=function(x,level=0.95)
{
  if(level>=1|level<0){return("the confidence level should be between 0 and 1!")}
  alpha=1-level
  n=length(x)
  t_a=t.quantile(df=n-1,1-alpha/2)
  barx=mean(x)
  s=sd(x)
  cat("The sample mean is", barx,", sample standard deviation is", s, ", and sample size is", n,"\n")
  cat("A", level*100, "% two-sided confidence interval for the population mean is [", barx-t_a*s/sqrt(n),",",bars+t_a*s/sqrt(n),"]","\n")
  t_a=t.quantile(df=n-1,1-alpha)
  cat("A", level*100, "% upper-confidence bound for the population mean is", barx+t_a*s/sqrt(n),"\n")
  cat("A", level*100, "% lower-confidence bound for the population mean is", barx-t_a*s/sqrt(n),"\n")
}

Tinterval.stat=function(n,barx,s,level=0.95)
{
  if(level>=1|level<0){return("the confidence level should be between 0 and 1!")}
  alpha=1-level
  t_a=t.quantile(df=n-1,1-alpha/2)
  cat("The sample mean is", barx,", sample standard deviation is", s, ", and sample size is", n,"\n")
  cat("A", level*100, "% two-sided confidence interval for the population mean is [", barx-t_a*s/sqrt(n),",",barx+t_a*s/sqrt(n),"]","\n")
  t_a=t.quantile(df=n-1,1-alpha)
  cat("A", level*100, "% upper-confidence bound for the population mean is", barx+t_a*s/sqrt(n),"\n")
  cat("A", level*100, "% lower-confidence bound for the population mean is", barx-t_a*s/sqrt(n),"\n")
}
