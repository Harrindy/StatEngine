Ztest=function(mu0,H1="two",alpha=0.05,sigma,
               data,n,barx)
{
  ########################################
  # Check input
  if(missing(mu0)==TRUE){return("Please provide the hypothesized value mu0 in the null hypothesis H0")}
  if(alpha>=1|alpha<=0){return("the significance level alpha should be between 0 and 1!")}
  if(missing(sigma)==TRUE){
    return("Please provide the value of sigma. If it is unknown, consider Ttest")
  }
  if(missing(data)==FALSE){
    n=length(data)
    barx=mean(data)
  }else if(missing(n)==TRUE | missing(barx)==TRUE)
  {
    return("please input the sample size/sample mean")
  }
  ########################################
  cat("The sample mean is", barx,"and sample size is", n,"\n","\n")
  z0=(barx-mu0)/(sigma/sqrt(n))
  if(H1=="two")
  {
  z_a=qnorm(1-alpha/2)
  cat("H1 is two-sided: mu does not equal to mu0=",mu0, ". The results are:","\n")
  if(abs(z0)<=z_a){
    cat("\n")
    cat("1. Test statistic z0 is", z0,", z_(alpha/2) is", z_a,
        ", Because |z0|<=z_(alpha/2), we fail to reject H0 at significance level", alpha,"\n")
  }else{
    cat("\n")
    cat("1. Test statistic z0 is", z0,", z_(alpha/2) is", z_a,
        ", Because |z0|>z_(alpha/2), we reject H0 at significance level", alpha,"\n")
  }
  Ln=barx-z_a*sigma/sqrt(n)
  Un=barx+z_a*sigma/sqrt(n)
  if(Ln<=mu0 & mu0<=Un){
    cat("\n")
    cat("2. A", (1-alpha)*100, "% two-sided confidence interval for the population mean is [", Ln,",",Un,"]","which contains the hypothesized value mu0=", mu0, ", we fail to reject H0 at significance level", alpha,"\n")
  }else{
    cat("\n")
    cat("2. A", (1-alpha)*100, "% two-sided confidence interval for the population mean is [", Ln,",",Un,"]","which does not contain the hypothesized value mu0=", mu0, ", we reject H0 at significance level", alpha,"\n")
  }
  pv=2*(1-pnorm(abs(z0)))
  if(pv< alpha){
    cat("\n")
    cat("3. The P-value is", pv ,"which is smaller than alpha=",alpha,", so we reject H0 at significance level", alpha,"\n")
  }else{
    cat("\n")
    cat("3. The P-value is", pv ,"which is not smaller than alpha=",alpha,", so we fail to reject H0 at significance level", alpha,"\n")
  }
  }else if(H1=="left")
  {
    z_a=qnorm(1-alpha)
    cat("H1 is one-sided: mu is less than mu0=",mu0, ". The results are:","\n")
    if(z0>= -z_a){
      cat("\n")
      cat("1. Test statistic z0 is", z0,", z_alpha is", z_a,
          ", Because -z_alpha<=z0, we fail to reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("1. Test statistic z0 is", z0,", z_alpha is", z_a,
          ", Because z0<-z_alpha, we reject H0 at significance level", alpha,"\n")
    }
    Un=barx+z_a*sigma/sqrt(n)
    if(mu0<=Un){
      cat("\n")
      cat("2. A", (1-alpha)*100, "% one-sided confidence interval for the population mean is (", -Inf,",",Un,"]","which contains the hypothesized value mu0=", mu0, ", we fail to reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("2. A", (1-alpha)*100, "% one-sided confidence interval for the population mean is (", -Inf,",",Un,"]","which does not contain the hypothesized value mu0=", mu0, ", we reject H0 at significance level", alpha,"\n")
    }
    pv=pnorm(z0)
    if(pv< alpha){
      cat("\n")
      cat("3. The P-value is", pv ,"which is smaller than alpha=",alpha,", so we reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("3. The P-value is", pv ,"which is not smaller than alpha=",alpha,", so we fail to reject H0 at significance level", alpha,"\n")
    }
  }else if(H1=="right")
  {
    z_a=qnorm(1-alpha)
    cat("H1 is one-sided: mu is more than mu0=",mu0, ". The results are:","\n")
    if(z0<= z_a){
      cat("\n")
      cat("1. Test statistic z0 is", z0,", z_alpha is", z_a,
          ", Because z0<= z_alpha, we fail to reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("1. Test statistic z0 is", z0,", z_alpha is", z_a,
          ", Because z0>z_alpha, we reject H0 at significance level", alpha,"\n")
    }
    Ln=barx-z_a*sigma/sqrt(n)
    if(mu0>=Ln){
      cat("\n")
      cat("2. A", (1-alpha)*100, "% one-sided confidence interval for the population mean is [", Ln,",",Inf,")","which contains the hypothesized value mu0=", mu0, ", we fail to reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("2. A", (1-alpha)*100, "% one-sided confidence interval for the population mean is [", Ln,",",Inf,")","which does not contain the hypothesized value mu0=", mu0, ", we reject H0 at significance level", alpha,"\n")
    }
    pv=1-pnorm(z0)
    if(pv< alpha){
      cat("\n")
      cat("3. The P-value is", pv ,"which is smaller than alpha=",alpha,", so we reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("3. The P-value is", pv ,"which is not smaller than alpha=",alpha,", so we fail to reject H0 at significance level", alpha,"\n")
    }
  }else{
    return("H1 should be two, left, or right")
  }
}

Ztest.power=function(H1="two",sigma,alpha=0.05,n,delta)
{
  ########################################
  # Check input
  if(ceiling(n)-n!=0 |n<=0){return("n should be a positive integer!")}
  if(alpha>=1|alpha<=0){return("the significance level alpha should be between 0 and 1!")}
  if(missing(sigma)==TRUE){
    return("Please provide the value of sigma. If it is unknown, consider Ttest")
  }
  ########################################
  if(H1=="two")
  {
    z_a=qnorm(1-alpha/2)
    cat("H1 is two-sided", "\n")
    if(missing(delta)!=TRUE){
      cat("\n")
      cat("The probability of the Type II error of this test at delta=mu1-mu0=",delta,"is",
          normal.prob(0,1,-z_a-delta*sqrt(n)/sigma,z_a-delta*sqrt(n)/sigma), "and the associated power is",
          1-normal.prob(0,1,-z_a-delta*sqrt(n)/sigma,z_a-delta*sqrt(n)/sigma))
    }
    par(mfrow=c(1,1))
    par(mar=c(4,4,4,.5))
    power=function(delta){1-normal.prob(0,1,-z_a-delta*sqrt(n)/sigma,z_a-delta*sqrt(n)/sigma)}
    delta=seq((-4-z_a)*sigma/sqrt(n),(4+z_a)*sigma/sqrt(n),length=100)
    pp=sapply(delta,power)
    plot(delta,pp,type="l",ylim=c(0,1),xlab=expression(delta==mu[1]-mu[0]),ylab="Power",main=paste("Power curve of the two-sided Z test at n=",ceiling(n),sep=""))
    lines(delta*2,rep(alpha,100),lty=3)
  }else if(H1=="left")
  {
    z_a=qnorm(1-alpha)
    cat("H1 is one-sided (<)","\n")
    if(missing(delta)!=TRUE){
      cat("\n")
      cat("The probability of the Type II error of this test at delta=mu1-mu0=",delta,"is",
          normal.prob(0,1,-z_a-delta*sqrt(n)/sigma,Inf), "and the associated power is",
          1-normal.prob(0,1,-z_a-delta*sqrt(n)/sigma,Inf))
    }
    par(mfrow=c(1,1))
    par(mar=c(4,4,4,.5))
    power=function(delta){1-normal.prob(0,1,-z_a-delta*sqrt(n)/sigma,Inf)}
    delta=seq((-4-z_a)*sigma/sqrt(n),0,length=100)
    pp=sapply(delta,power)
    plot(delta,pp,type="l",ylim=c(0,1),xlab=expression(delta==mu[1]-mu[0]),ylab="Power",main=paste("Power curve of the One-sided (<) Z test at n=",ceiling(n),sep=""))
    lines(delta*2,rep(alpha,100),lty=3)
  }else if(H1=="right")
  {
    z_a=qnorm(1-alpha)
    cat("H1 is one-sided (>)","\n")
    if(missing(delta)!=TRUE){
      cat("\n")
      cat("The probability of the Type II error of this test at delta=mu1-mu0=",delta,"is",
          normal.prob(0,1,-Inf,z_a-delta*sqrt(n)/sigma), "and the associated power is",
          1-normal.prob(0,1,-Inf,z_a-delta*sqrt(n)/sigma))
    }
    par(mfrow=c(1,1))
    par(mar=c(4,4,4,.5))
    power=function(delta){1-normal.prob(0,1,-Inf,z_a-delta*sqrt(n)/sigma)}
    delta=seq(0,(4+z_a)*sigma/sqrt(n),length=100)
    pp=sapply(delta,power)
    plot(delta,pp,type="l",ylim=c(0,1),xlab=expression(delta==mu[1]-mu[0]),ylab="Power",main=paste("Power curve of the one-sided (>) Z test at n=",ceiling(n),sep=""))
    lines(delta*2,rep(alpha,100),lty=3)
  }else{
    return("H1 should be two, right, or right.")
  }
}

Ztest.sample.size=function(H1="two",alpha,beta,sigma,delta)
{
  if(alpha>=1|alpha<=0){return("the significance level alpha should be between 0 and 1!")}
  if(beta>=1|beta<=0){return("the require Type II error should be of a probability beta between 0 and 1!")}
  if(missing(sigma)==TRUE){
    return("Please provide the value of sigma. If it is unknown, consider Ttest")
  }
  if(missing(delta)==TRUE){
    return("Please provide the value of delta")
  }
  if(H1=="left" & delta>0){return("H1 only holds when delta<0")}
  if(H1=="right" & delta<0){return("H1 only holds when delta>0")}
  if(H1=="two"){
    z_a=qnorm(1-alpha/2)
    pp=function(n){
      normal.prob(0,1,-z_a-delta*sqrt(n)/sigma,z_a-delta*sqrt(n)/sigma)
    }
    targ=function(n){return((pp(n)-beta)^2)}
    n=optimize(targ,c(1,4*ceiling(((4-z_a)*sigma/delta)^2)))$minimum
    cat("At significance level alpha=",alpha,", we need at least n=",ceiling(n),"to make the power of this test at delta=",delta,"be at least", 1-beta)
  }else if(H1=="left"){
    z_a=qnorm(1-alpha)
    z_beta=qnorm(beta)
    n=((z_a-z_beta)*sigma/delta)^2
    cat("At significance level alpha=",alpha,", we need at least n=",ceiling(n),"to make the power of this test at delta=",delta,"be at least", 1-beta)
  }else if(H1=="right"){
    z_a=qnorm(1-alpha)
    z_beta=qnorm(1-beta)
    n=((-z_a-z_beta)*sigma/delta)^2
    cat("At significance level alpha=",alpha,", we need at least n=",ceiling(n),"to make the power of this test at delta=",delta,"be at least", 1-beta)
  }else{
    return("H1 should be two, right, or right.")
  }
}



Ttest=function(mu0,H1="two",alpha=0.05,
               data,n,barx,s)
{
  ########################################
  # Check input
  if(missing(mu0)==TRUE){return("Please provide the hypothesized value mu0 in the null hypothesis H0")}
  if(alpha>=1|alpha<=0){return("the significance level alpha should be between 0 and 1!")}
  if(missing(data)==FALSE){
    n=length(data)
    barx=mean(data)
    s=sd(data)
  }else if(missing(n)==TRUE | missing(barx)==TRUE)
  {
    return("please input the sample size/sample mean")
  }
  ########################################
  cat("The sample mean is", barx,"sample standard deviation is", s,", and sample size is", n,"\n","\n")
  t0=(barx-mu0)/(s/sqrt(n))
  if(H1=="two")
  {
    t_a=t.quantile(n-1,1-alpha/2)
    cat("H1 is two-sided: mu does not equal to mu0=",mu0, ". The results are:","\n")
    if(abs(t0)<= t_a){
      cat("\n")
      cat("1. Test statistic t0 is", t0,", t_(n-1,alpha/2) is", t_a,
          ", Because |t0|<=t_(n-1,alpha/2), we fail to reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("1. Test statistic t0 is", t0,", t_(n-1,alpha/2) is", t_a,
          ", Because |t0|>t_(n-1,alpha/2), we reject H0 at significance level", alpha,"\n")
    }
    Ln=barx-t_a*s/sqrt(n)
    Un=barx+t_a*s/sqrt(n)
    if(Ln<=mu0 & mu0<=Un){
      cat("\n")
      cat("2. A", (1-alpha)*100, "% two-sided confidence interval for the population mean is [", Ln,",",Un,"]","which contains the hypothesized value mu0=", mu0, ", we fail to reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("2. A", (1-alpha)*100, "% two-sided confidence interval for the population mean is [", Ln,",",Un,"]","which does not contain the hypothesized value mu0=", mu0, ", we reject H0 at significance level", alpha,"\n")
    }
    pv=2*(1-pt(abs(t0),df=n-1))
    if(pv< alpha){
      cat("\n")
      cat("3. The P-value is", pv ,"which is smaller than alpha=",alpha,", so we reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("3. The P-value is", pv ,"which is not smaller than alpha=",alpha,", so we fail to reject H0 at significance level", alpha,"\n")
    }
  }else if(H1=="left")
  {
    t_a=t.quantile(n-1,1-alpha)
    cat("H1 is one-sided: mu is less than mu0=",mu0, ". The results are:","\n")
    if(t0>= -t_a){
      cat("\n")
      cat("1. Test statistic t0 is", t0,", t_(n-1,alpha) is", t_a,
          ", Because -t_(n-1,alpha)<=t0, we fail to reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("1. Test statistic t0 is", t0,", t_(n-1,alpha) is", t_a,
          ", Because t0<-t_(n-1,alpha), we reject H0 at significance level", alpha,"\n")
    }
    Un=barx+t_a*s/sqrt(n)
    if(mu0<=Un){
      cat("\n")
      cat("2. A", (1-alpha)*100, "% one-sided confidence interval for the population mean is (", -Inf,",",Un,"]","which contains the hypothesized value mu0=", mu0, ", we fail to reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("2. A", (1-alpha)*100, "% one-sided confidence interval for the population mean is (", -Inf,",",Un,"]","which does not contain the hypothesized value mu0=", mu0, ", we reject H0 at significance level", alpha,"\n")
    }
    pv=pt(t0,df=n-1)
    if(pv< alpha){
      cat("\n")
      cat("3. The P-value is", pv ,"which is smaller than alpha=",alpha,", so we reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("3. The P-value is", pv ,"which is not smaller than alpha=",alpha,", so we fail to reject H0 at significance level", alpha,"\n")
    }
  }else if(H1=="right")
  {
    t_a=t.quantile(n-1,1-alpha)
    cat("H1 is one-sided: mu is more than mu0=",mu0, ". The results are:","\n")
    if(t0<= t_a){
      cat("\n")
      cat("1. Test statistic t0 is", t0,", t_(n-1,alpha) is", t_a,
          ", Because t0<= t_(n-1,alpha), we fail to reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("1. Test statistic t0 is", t0,", t_(n-1,alpha) is", t_a,
          ", Because t0>t_(n-1,alpha), we reject H0 at significance level", alpha,"\n")
    }
    Ln=barx-t_a*s/sqrt(n)
    if(mu0>=Ln){
      cat("\n")
      cat("2. A", (1-alpha)*100, "% one-sided confidence interval for the population mean is [", Ln,",",Inf,")","which contains the hypothesized value mu0=", mu0, ", we fail to reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("2. A", (1-alpha)*100, "% one-sided confidence interval for the population mean is [", Ln,",",Inf,")","which does not contain the hypothesized value mu0=", mu0, ", we reject H0 at significance level", alpha,"\n")
    }
    pv=1-pt(t0,df=n-1)
    if(pv< alpha){
      cat("\n")
      cat("3. The P-value is", pv ,"which is smaller than alpha=",alpha,", so we reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("3. The P-value is", pv ,"which is not smaller than alpha=",alpha,", so we fail to reject H0 at significance level", alpha,"\n")
    }
  }else{
    return("H1 should be two, left, or right")
  }
}
