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
    cat("2. A", (1-alpha)*100, "% two-sided confidence interval for the population mean is [", Ln,",",Un,"]","which contains the hypothesized value mu0=", mu0, ", so we fail to reject H0 at significance level", alpha,"\n")
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
      cat("2. A", (1-alpha)*100, "% one-sided confidence interval for the population mean is (", -Inf,",",Un,"]","which contains the hypothesized value mu0=", mu0, ", so we fail to reject H0 at significance level", alpha,"\n")
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
      cat("2. A", (1-alpha)*100, "% one-sided confidence interval for the population mean is [", Ln,",",Inf,")","which contains the hypothesized value mu0=", mu0, ", so we fail to reject H0 at significance level", alpha,"\n")
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
    cat("H1 is left-tailed (<)","\n")
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
    cat("H1 is right-tailed (>)","\n")
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
    plot(delta,pp,type="l",ylim=c(0,1),xlab=expression(delta==mu[1]-mu[0]),ylab="Power",main=paste("Power curve of the right-tailed (>) Z test at n=",ceiling(n),sep=""))
    lines(delta*2,rep(alpha,100),lty=3)
  }else{
    return("H1 should be two, left, or right.")
  }
}

sample.size.Ztest=function(H1="two",alpha,beta,sigma,delta)
{
  if(alpha>=1|alpha<=0){return("The significance level alpha should be between 0 and 1!")}
  if(beta>=1|beta<=alpha){return(cat("The specified Type II error should be of a probability beta between alpha=",alpha, "and 1!"))}
  if(missing(sigma)==TRUE){
    return("Please provide the value of sigma. If it is unknown, consider Ttest")
  }
  if(missing(delta)==TRUE){
    return("Please provide the value of delta")
  }
  if(delta==0){return("delta should not be 0")}
  if(H1=="left" & delta>=0){return("H1 only holds when delta<0")}
  if(H1=="right" & delta<=0){return("H1 only holds when delta>0")}
  if(H1=="two"){
    z_a=qnorm(1-alpha/2)
    pp=function(n){
      normal.prob(0,1,-z_a-delta*sqrt(n)/sigma,z_a-delta*sqrt(n)/sigma)
    }
    targ=function(n){return(pp(n)-beta)}
    n=ceiling(uniroot(targ,c(1,4*ceiling(((4-z_a)*sigma/delta)^2)))$root)
    cat("At significance level alpha=",alpha,", we need at least n=",n,"to achieve a power >=",1-beta, "of this test at delta=",delta,"\n")
    cat("When n=",n-1,", the power is", 1-pp(n-1),"\n")
    cat("When n=",n,", the power is", 1-pp(n),"\n")
    cat("When n=",n+1,", the power is", 1-pp(n+1),"\n")
  }else if(H1=="left"){
    z_a=qnorm(1-alpha)
    z_beta=qnorm(beta)
    pp=function(n){
      normal.prob(0,1,-z_a-delta*sqrt(n)/sigma,Inf)
    }
    n=ceiling(((z_a-z_beta)*sigma/delta)^2)
    cat("At significance level alpha=",alpha,", we need at least n=",n,"to achieve a power >=",1-beta, "of this test at delta=",delta,"\n")
    cat("When n=",n-1,", the power is", 1-pp(n-1),"\n")
    cat("When n=",n,", the power is", 1-pp(n),"\n")
    cat("When n=",n+1,", the power is", 1-pp(n+1),"\n")
    }else if(H1=="right"){
    z_a=qnorm(1-alpha)
    z_beta=qnorm(1-beta)
    pp=function(n){
      normal.prob(0,1,-Inf,z_a-delta*sqrt(n)/sigma)
    }
    n=ceiling(((-z_a-z_beta)*sigma/delta)^2)
    cat("At significance level alpha=",alpha,", we need at least n=",n,"to achieve a power >=",1-beta, "of this test at delta=",delta,"\n")
    cat("When n=",n-1,", the power is", 1-pp(n-1),"\n")
    cat("When n=",n,", the power is", 1-pp(n),"\n")
    cat("When n=",n+1,", the power is", 1-pp(n+1),"\n")
  }else{
    return("H1 should be two, left, or right.")
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
      cat("2. A", (1-alpha)*100, "% two-sided confidence interval for the population mean is [", Ln,",",Un,"]","which contains the hypothesized value mu0=", mu0, ", so we fail to reject H0 at significance level", alpha,"\n")
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
      cat("2. A", (1-alpha)*100, "% one-sided confidence interval for the population mean is (", -Inf,",",Un,"]","which contains the hypothesized value mu0=", mu0, ", so we fail to reject H0 at significance level", alpha,"\n")
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
      cat("2. A", (1-alpha)*100, "% one-sided confidence interval for the population mean is [", Ln,",",Inf,")","which contains the hypothesized value mu0=", mu0, ", so we fail to reject H0 at significance level", alpha,"\n")
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


Ttest.power=function(H1="two",est.sigma,alpha=0.05,n,delta)
{
  ########################################
  # Check input
  if(ceiling(n)-n!=0 |n<=0){return("n should be a positive integer!")}
  if(alpha>=1|alpha<=0){return("the significance level alpha should be between 0 and 1!")}
  if(missing(est.sigma)==TRUE){
    return("Please provide an estimate of sigma.")
  }
  ########################################
  if(H1=="two")
  {
    t_a=t.quantile(n-1,1-alpha/2)
    beta=pt(t_a,df=n-1,ncp=delta*sqrt(n)/est.sigma)-pt(-t_a,df=n-1,ncp=delta*sqrt(n)/est.sigma)
    cat("H1 is two-sided", "\n")
    if(missing(delta)!=TRUE){
      cat("\n")
      cat("The probability of the Type II error of this test at delta=mu1-mu0=",delta,"is",
          beta, "and the associated power is",
          1-beta)
    }
    par(mfrow=c(1,1))
    par(mar=c(4,4,4,.5))
    power=function(delta){
      beta=pt(t_a,df=n-1,ncp=abs(delta)*sqrt(n)/est.sigma)-pt(-t_a,df=n-1,ncp=abs(delta)*sqrt(n)/est.sigma)
      return(1-beta)
    }
    targ=function(delta){
      res=(power(delta)-0.9995)
      return(res)
    }
    a=uniroot(targ,c(0,100))$root
    delta=seq(0,a,length=100)
    pp=sapply(delta,power)
    plot(delta,pp,type="l",ylim=c(0,1),xlim=c(-a,a),xlab=expression(delta==mu[1]-mu[0]),ylab="Power",main=paste("Power curve of the two-sided T test at n=",ceiling(n),sep=""))
    lines(-delta,pp)
    lines(seq(-2*a,2*a,length=100),rep(alpha,100),lty=3)
  }else if(H1=="left")
  {
    t_a=t.quantile(n-1,1-alpha)
    cat("H1 is left-tailed (<)","\n")
    beta=1-pt(-t_a,df=n-1,ncp=delta*sqrt(n)/est.sigma)
    if(missing(delta)!=TRUE){
      cat("\n")
      cat("The probability of the Type II error of this test at delta=mu1-mu0=",delta,"is",
          beta, "and the associated power is",
          1-beta)
    }
    par(mfrow=c(1,1))
    par(mar=c(4,4,4,.5))
    power=function(delta){pt(-t_a,df=n-1,ncp=delta*sqrt(n)/est.sigma)}
    targ=function(delta){
      res=(power(delta)-0.9995)
      return(res)
    }
    a=uniroot(targ,c(-10,0))$root
    delta=seq(a,0,length=100)
    pp=sapply(delta,power)
    plot(delta,pp,type="l",ylim=c(0,1),xlab=expression(delta==mu[1]-mu[0]),ylab="Power",main=paste("Power curve of the left-tailed (<) T test at n=",ceiling(n),sep=""))
    lines(delta*2,rep(alpha,100),lty=3)
  }else if(H1=="right")
  {
    t_a=t.quantile(n-1,1-alpha)
    cat("H1 is left-tailed (<)","\n")
    beta=pt(t_a,df=n-1,ncp=delta*sqrt(n)/est.sigma)
    if(missing(delta)!=TRUE){
      cat("\n")
      cat("The probability of the Type II error of this test at delta=mu1-mu0=",delta,"is",
          beta, "and the associated power is",
          1-beta)
    }
    par(mfrow=c(1,1))
    par(mar=c(4,4,4,.5))
    power=function(delta){1-pt(t_a,df=n-1,ncp=delta*sqrt(n)/est.sigma)}
    targ=function(delta){
      res=(power(delta)-0.9995)
      return(res)
    }
    a=uniroot(targ,c(0,10))$root
    delta=seq(a,0,length=100)
    pp=sapply(delta,power)
    plot(delta,pp,type="l",ylim=c(0,1),xlab=expression(delta==mu[1]-mu[0]),ylab="Power",main=paste("Power curve of the right-tailed (>) T test at n=",ceiling(n),sep=""))
    lines(delta*2,rep(alpha,100),lty=3)
  }else{
    return("H1 should be two, left, or right.")
  }
}


sample.size.Ttest=function(H1="two",est.sigma,beta,delta,alpha=0.05)
{
  if(alpha>=1|alpha<=0){return("The significance level alpha should be between 0 and 1!")}
  if(beta>=1|beta<=alpha){return(cat("The specified Type II error should be of a probability beta between alpha=",alpha, "and 1!"))}
  if(missing(est.sigma)==TRUE){
    return("Please provide an estimate of sigma.")
  }
  if(missing(delta)==TRUE){
    return("Please provide the value of delta")
  }
  if(delta==0){return("delta should not be 0")}
  if(H1=="left" & delta>=0){return("H1 only holds when delta<0")}
  if(H1=="right" & delta<=0){return("H1 only holds when delta>0")}
  if(H1=="two"){
    pp=function(n){
      t_a=t.quantile(n-1,1-alpha/2)
      return(pt(t_a,df=n-1,ncp=abs(delta)*sqrt(n)/est.sigma)-pt(-t_a,df=n-1,ncp=abs(delta)*sqrt(n)/est.sigma))
    }
    targ=function(n){return(pp(n)-beta)}
    z_a=normal.quantile(0,1,alpha/2)
    n=ceiling(uniroot(targ, c(2,4*ceiling(((4-z_a)*est.sigma/delta)^2)) )$root)
    cat("At significance level alpha=",alpha,", we need at least n=",n,"to achieve a power >=",1-beta, "of this test at delta=",delta,"\n")
    cat("When n=",n-1,", the power is", 1-pp(n-1),"\n")
    cat("When n=",n,", the power is", 1-pp(n),"\n")
    cat("When n=",n+1,", the power is", 1-pp(n+1),"\n")
  }else if(H1=="left"){
    pp=function(n){
      t_a=t.quantile(n-1,1-alpha/2)
      return(1-pt(-t_a,df=n-1,ncp=delta*sqrt(n)/est.sigma))
    }
    z_a=normal.quantile(0,1,alpha/2)
    targ=function(n){return(pp(n)-beta)}
    n=ceiling(uniroot(targ, c(2,4*ceiling(((4-z_a)*est.sigma/delta)^2)) )$root)
    cat("At significance level alpha=",alpha,", we need at least n=",n,"to achieve a power >=",1-beta, "of this test at delta=",delta,"\n")
    cat("When n=",n-1,", the power is", 1-pp(n-1),"\n")
    cat("When n=",n,", the power is", 1-pp(n),"\n")
    cat("When n=",n+1,", the power is", 1-pp(n+1),"\n")
  }else if(H1=="right"){
    pp=function(n){
      t_a=t.quantile(n-1,1-alpha/2)
      return(pt(t_a,df=n-1,ncp=delta*sqrt(n)/est.sigma))
    }
    targ=function(n){return(pp(n)-beta)}
    n=ceiling(uniroot(targ, c(2,4*ceiling(((4-z_a)*est.sigma/delta)^2)) )$root)
    cat("At significance level alpha=",alpha,", we need at least n=",n,"to achieve a power >=",1-beta, "of this test at delta=",delta,"\n")
    cat("When n=",n-1,", the power is", 1-pp(n-1),"\n")
    cat("When n=",n,", the power is", 1-pp(n),"\n")
    cat("When n=",n+1,", the power is", 1-pp(n+1),"\n")
  }else{
    return("H1 should be two, left, or right.")
  }
}
