Zinterval.data=function(data,sigma,level=0.95)
{
  if(level>=1|level<0){return("the confidence level should be between 0 and 1!")}
  alpha=1-level
  z_a=qnorm(1-alpha/2)
  n=length(data)
  barx=mean(data)
  cat("The sample mean is", mean(data),"and sample size is", n,"\n")
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
  cat("In the estimation of mu, the samllest sample size to control the margin of error E <=", E, "at",level*100, "% confidence level is", n, "\n")
}

AZinterval.data=function(data,level=0.95)
{
  s=sd(data)
  if(level>=1|level<0){return("the confidence level should be between 0 and 1!")}
  alpha=1-level
  z_a=qnorm(1-alpha/2)
  n=length(data)
  barx=mean(data)
  cat("The sample mean is", mean(data),", sample standard deviation is", s, ", and sample size is", n,"\n")
  cat("A large-sample confidence interval for the population mean with confidence level of approximately", level*100, "% is [", barx-z_a*s/sqrt(n),",",barx+z_a*s/sqrt(n),"]","\n")
  z_a=qnorm(1-alpha)
  cat("A large-sample upper-confidence bound for the population mean with confidence level of approximately", level*100, "% is", barx+z_a*s/sqrt(n),"\n")
  cat("A large-sample lower-confidence bound for the population mean with confidence level of approximately", level*100, "% is", barx-z_a*s/sqrt(n),"\n")
}

AZinterval.stat=function(n,barx,s,level=0.95)
{
  if(level>=1|level<0){return("the confidence level should be between 0 and 1!")}
  alpha=1-level
  z_a=qnorm(1-alpha/2)
  cat("The sample mean is", barx,", sample standard deviation is", s, ", and sample size is", n,"\n")
  cat("A large-sample confidence interval for the population mean with confidence level of approximately", level*100, "% is [", barx-z_a*s/sqrt(n),",",barx+z_a*s/sqrt(n),"]","\n")
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


Tinterval.data=function(data,level=0.95)
{
  if(level>=1|level<0){return("the confidence level should be between 0 and 1!")}
  alpha=1-level
  n=length(data)
  t_a=t.quantile(df=n-1,1-alpha/2)
  barx=mean(data)
  s=sd(data)
  cat("The sample mean is", barx,", sample standard deviation is", s, ", and sample size is", n,"\n")
  cat("A", level*100, "% two-sided confidence interval for the population mean is [", barx-t_a*s/sqrt(n),",",barx+t_a*s/sqrt(n),"]","\n")
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

Chi2.quantile=function(df,q)
{
  if(df<=0){return("df must be positive")}
  if(q<=0|q>=1){return("q must be between 0 and 1")}
  return(qchisq(q,df))
}

Chi2interval.data=function(data,level=0.95)
{
  if(level>=1|level<0){return("the confidence level should be between 0 and 1!")}
  alpha=1-level
  n=length(data)
  chi_a1=Chi2.quantile(df=n-1,1-alpha/2)
  chi_a2=Chi2.quantile(df=n-1,alpha/2)
  s2=var(data)
  cat("The sample standard variance is", s2, "and sample size is", n,"\n")
  cat("A", level*100, "% two-sided confidence interval for the population variance is [", (n-1)*s2/chi_a1,",",(n-1)*s2/chi_a2,"]","\n")
  chi_a1=Chi2.quantile(df=n-1,1-alpha)
  chi_a2=Chi2.quantile(df=n-1,alpha)
  cat("A", level*100, "% upper-confidence bound for the population variance is", (n-1)*s2/chi_a2,"\n")
  cat("A", level*100, "% lower-confidence bound for the population variance is", (n-1)*s2/chi_a1,"\n","\n")

  chi_a1=Chi2.quantile(df=n-1,1-alpha/2)
  chi_a2=Chi2.quantile(df=n-1,alpha/2)
  cat("The sample standard deviation is", sqrt(s2), "and sample size is", n,"\n")
  cat("A", level*100, "% two-sided confidence interval for the population standard deviation is [", sqrt((n-1)*s2/chi_a1),",",sqrt((n-1)*s2/chi_a2),"]","\n")
  chi_a1=Chi2.quantile(df=n-1,1-alpha)
  chi_a2=Chi2.quantile(df=n-1,alpha)
  cat("A", level*100, "% upper-confidence bound for the population standard deviation is", sqrt((n-1)*s2/chi_a2),"\n")
  cat("A", level*100, "% lower-confidence bound for the population standard deviation is", sqrt((n-1)*s2/chi_a1),"\n")
}

Chi2interval.stat=function(n,s,level=0.95)
{
  if(level>=1|level<0){return("the confidence level should be between 0 and 1!")}
  alpha=1-level
  chi_a1=Chi2.quantile(df=n-1,1-alpha/2)
  chi_a2=Chi2.quantile(df=n-1,alpha/2)
  s2=s^2
  cat("The sample standard variance is", s2, "and sample size is", n,"\n")
  cat("A", level*100, "% two-sided confidence interval for the population variance is [", (n-1)*s2/chi_a1,",",(n-1)*s2/chi_a2,"]","\n")
  chi_a1=Chi2.quantile(df=n-1,1-alpha)
  chi_a2=Chi2.quantile(df=n-1,alpha)
  cat("A", level*100, "% upper-confidence bound for the population variance is", (n-1)*s2/chi_a2,"\n")
  cat("A", level*100, "% lower-confidence bound for the population variance is", (n-1)*s2/chi_a1,"\n","\n")

  chi_a1=Chi2.quantile(df=n-1,1-alpha/2)
  chi_a2=Chi2.quantile(df=n-1,alpha/2)
  cat("The sample standard deviation is", sqrt(s2), "and sample size is", n,"\n")
  cat("A", level*100, "% two-sided confidence interval for the population standard deviation is [", sqrt((n-1)*s2/chi_a1),",",sqrt((n-1)*s2/chi_a2),"]","\n")
  chi_a1=Chi2.quantile(df=n-1,1-alpha)
  chi_a2=Chi2.quantile(df=n-1,alpha)
  cat("A", level*100, "% upper-confidence bound for the population standard deviation is", sqrt((n-1)*s2/chi_a2),"\n")
  cat("A", level*100, "% lower-confidence bound for the population standard deviation is", sqrt((n-1)*s2/chi_a1),"\n")
}

Propinterval=function(n,x,level=0.95)
{
  if(level>=1|level<0){return("the confidence level should be between 0 and 1!")}
  p=x/n
  alpha=1-level
  z_a=qnorm(1-alpha/2)
  cat("The sample proportion is", p,"and sample size is", n,"\n")
  cat("A large-sample confidence interval for the population proportion with confidence level of approximately", level*100, "% is [", p-z_a*sqrt(p*(1-p)/n),",",p+z_a*sqrt(p*(1-p)/n),"]","\n")
  z_a=qnorm(1-alpha)
  cat("A large-sample upper-confidence bound for the population proportion with confidence level of approximately", level*100, "% is", p+z_a*sqrt(p*(1-p)/n),"\n")
  cat("A large-sample lower-confidence bound for the population proportion with confidence level of approximately", level*100, "% is", p-z_a*sqrt(p*(1-p)/n),"\n")
}

sample.size.Propinterval=function(E,ini.p=0.5,level=0.95)
{
  if(level>=1|level<0){return("the confidence level should be between 0 and 1!")}
  if(ini.p>=1|level<0){return("the initial estimate should be between 0 and 1!")}
  n=ceiling((qnorm(1-(1-level)/2)/E)^2*ini.p*(1-ini.p))
  cat("In the estimation of p, the samllest sample size to control the margin of error E <=", E, "at",level*100, "% confidence level is", n, "\n")
}

Predinterval.data=function(data,level=0.95)
{
  if(level>=1|level<0){return("the confidence level should be between 0 and 1!")}
  alpha=1-level
  n=length(data)
  t_a=t.quantile(df=n-1,1-alpha/2)
  barx=mean(data)
  s=sd(data)
  cat("The sample mean is", barx,", sample standard deviation is", s, ", and sample size is", n,"\n")
  cat("A", level*100, "% two-sided predict interval for the next observation is [", barx-t_a*s*sqrt(1+1/n),",",barx+t_a*s*sqrt(1+1/n),"]","\n")
  t_a=t.quantile(df=n-1,1-alpha)
  cat("A", level*100, "% upper-confidence bound for the next observation is", barx+t_a*s*sqrt(1+1/n),"\n")
  cat("A", level*100, "% lower-confidence bound for the next observation is", barx-t_a*s*sqrt(1+1/n),"\n")
}

Predinterval.stat=function(n,barx,s,level=0.95)
{
  if(level>=1|level<0){return("the confidence level should be between 0 and 1!")}
  alpha=1-level
  t_a=t.quantile(df=n-1,1-alpha/2)
  cat("The sample mean is", barx,", sample standard deviation is", s, ", and sample size is", n,"\n")
  cat("A", level*100, "% two-sided predict interval for the next observation is [", barx-t_a*s*sqrt(1+1/n),",",barx+t_a*s*sqrt(1+1/n),"]","\n")
  t_a=t.quantile(df=n-1,1-alpha)
  cat("A", level*100, "% upper-confidence bound for the next observation is", barx+t_a*s*sqrt(1+1/n),"\n")
  cat("A", level*100, "% lower-confidence bound for the next observation is", barx-t_a*s*sqrt(1+1/n),"\n")
}
