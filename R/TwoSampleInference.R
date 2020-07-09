twosample.Zinterval=function(level,sigma1,sigma2,sample1,sample2,barx1,barx2,n1,n2)
{
  if(level>=1|level<=0){return("the confidence level should be between 0 and 1!")}
  if(missing(sigma1)==TRUE|missing(sigma2)==TRUE){
    return("Please provide the value of sigma1 and sigma2.
           If they are unknown, consider twosample.Tinterval")
  }
  if(missing(sample1)==FALSE){
    n1=length(sample1)
    barx1=mean(sample1)
  }
  if(missing(sample2)==FALSE){
    n2=length(sample2)
    barx2=mean(sample2)
  }
  alpha=1-level
  z_a=qnorm(1-alpha/2)
  s=sqrt(sigma1^2/n1+sigma2^2/n2)
  cat("A", level*100, "% two-sided confidence interval for the difference in population means is [", barx1-barx2-z_a*s,",",barx1-barx2+z_a*s,"]","\n")
  z_a=qnorm(1-alpha)
  cat("A", level*100, "% upper-confidence bound for the population mean is", barx1-barx2+z_a*s,"\n")
  cat("A", level*100, "% lower-confidence bound for the population mean is", barx1-barx2-z_a*s,"\n")
}

twosample.Tinterval=function(level,pooled="no",sample1,sample2,barx1,barx2,n1,n2,s1,s2)
{
  if(level>=1|level<=0){return("the confidence level should be between 0 and 1!")}
  if(missing(sample1)==FALSE){
    n1=length(sample1)
    barx1=mean(sample1)
    s1=sd(sample1)
  }
  if(missing(sample2)==FALSE){
    n2=length(sample2)
    barx2=mean(sample2)
    s2=sd(sample2)
  }
  if(pooled=="yes"){
    v=n1+n2-2
    s=sqrt((n1-1)*s1^2+(n2-1)*s2^2)/sqrt(n1+n2-2)*sqrt(1/n1+1/n2)
  }else if(pooled=="no"){
    v=(s1^2/n1+s2^2/n2)^2/((s1^2/n1)^2/(n1-1)+(s2^2/n2)^2/(n2-1))
    s=sqrt(s1^2/n1+s2^2/n2)
  }
  alpha=1-level
  z_a=qt(1-alpha/2,df=v)
  cat("A", level*100, "% two-sided confidence interval for the difference in population means is [", barx1-barx2-z_a*s,",",barx1-barx2+z_a*s,"]","\n")
  z_a=qt(1-alpha,df=v)
  cat("A", level*100, "% upper-confidence bound for the population mean is", barx1-barx2+z_a*s,"\n")
  cat("A", level*100, "% lower-confidence bound for the population mean is", barx1-barx2-z_a*s,"\n")
}

twosample.Ztest=function(Delta0=0,H1="two",alpha,sigma1,sigma2,sample1,sample2,barx1,barx2,n1,n2)
{
  if(missing(Delta0)==TRUE){return("Please provide the hypothesized value Delta0 in the null hypothesis H0")}
  if(alpha>=1|alpha<=0){return("the significance level alpha should be between 0 and 1!")}
  if(missing(sigma1)==TRUE|missing(sigma2)==TRUE){
    return("Please provide the value of sigma1 and sigma2.
           If they are unknown, consider twosample.Tinterval")
  }
  if(missing(sample1)==FALSE){
    n1=length(sample1)
    barx1=mean(sample1)
  }
  if(missing(sample2)==FALSE){
    n2=length(sample2)
    barx2=mean(sample2)
  }
  z0=(barx1-barx2-Delta0)/sqrt(sigma1^2/n1+sigma2^2/n2)

  if(H1=="two"){
    z_a=qnorm(1-alpha/2)
    cat("H1 is two-tailed. The results are:","\n")
    if(abs(z0)<=z_a){
      cat("\n")
      cat("1. Test statistic z0 is", z0,", z_(alpha/2) is", z_a,
          ". Because |z0|<=z_(alpha/2), we fail to reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("1. Test statistic z0 is", z0,", z_(alpha/2) is", z_a,
          ". Because |z0|>z_(alpha/2), so we reject H0 at significance level", alpha,"\n")
    }
    Ln=barx1-barx2-z_a*sqrt(sigma1^2/n1+sigma2^2/n2)
    Un=barx1-barx2+z_a*sqrt(sigma1^2/n1+sigma2^2/n2)
    if(Ln<=Delta0 & Delta0<=Un){
      cat("\n")
      cat("2. A", (1-alpha)*100, "% two-tailed confidence interval for the population mean is [", Ln,",",Un,"]","which contains the hypothesized value Delta0=", Delta0, ", so we fail to reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("2. A", (1-alpha)*100, "% two-tailed confidence interval for the population mean is [", Ln,",",Un,"]","which does not contain the hypothesized value Delta0=", Delta0, ", so we reject H0 at significance level", alpha,"\n")
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
    cat("H1 is left-tailed. The results are:","\n")
    if(z0>= -z_a){
      cat("\n")
      cat("1. Test statistic z0 is", z0,", z_alpha is", z_a,
          ". Because -z_alpha<=z0, we fail to reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("1. Test statistic z0 is", z0,", z_alpha is", z_a,
          ". Because z0<-z_alpha, so we reject H0 at significance level", alpha,"\n")
    }
    Un=barx1-barx2+z_a*sqrt(sigma1^2/n1+sigma2^2/n2)
    if(Delta0<=Un){
      cat("\n")
      cat("2. A", (1-alpha)*100, "% one-sided confidence interval for the population mean is (", -Inf,",",Un,"]","which contains the hypothesized value Delta0=", Delta0, ", so we fail to reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("2. A", (1-alpha)*100, "% one-sided confidence interval for the population mean is (", -Inf,",",Un,"]","which does not contain the hypothesized value Delta0=", Delta0, ", so we reject H0 at significance level", alpha,"\n")
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
    cat("H1 is right-tailed. The results are:","\n")
    if(z0<= z_a){
      cat("\n")
      cat("1. Test statistic z0 is", z0,", z_alpha is", z_a,
          ". Because z0<= z_alpha, we fail to reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("1. Test statistic z0 is", z0,", z_alpha is", z_a,
          ". Because z0>z_alpha, so we reject H0 at significance level", alpha,"\n")
    }
    Ln=barx1-barx2-z_a*sqrt(sigma1^2/n1+sigma2^2/n2)
    if(Delta0>=Ln){
      cat("\n")
      cat("2. A", (1-alpha)*100, "% one-sided confidence interval for the population mean is [", Ln,",",Inf,")","which contains the hypothesized value Delta0=", Delta0, ", so we fail to reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("2. A", (1-alpha)*100, "% one-sided confidence interval for the population mean is [", Ln,",",Inf,")","which does not contain the hypothesized value Delta0=", Delta0, ", so we reject H0 at significance level", alpha,"\n")
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


twosample.Ttest=function(Delta0=0,H1="two",alpha,pooled="no",sample1,sample2,barx1,n1,s1,barx2,n2,s2)
{
  if(missing(Delta0)==TRUE){return("Please provide the hypothesized value Delta0 in the null hypothesis H0")}
  if(alpha>=1|alpha<=0){return("the significance level alpha should be between 0 and 1!")}
  if(missing(sample1)==FALSE){
    n1=length(sample1)
    barx1=mean(sample1)
    s1=sd(sample1)
  }
  if(missing(sample2)==FALSE){
    n2=length(sample2)
    barx2=mean(sample2)
    s2=sd(sample2)
  }
  if(pooled=="yes"){
    v=n1+n2-2
    s=sqrt((n1-1)*s1^2+(n2-1)*s2^2)/sqrt(n1+n2-2)*sqrt(1/n1+1/n2)
  }else if(pooled=="no"){
    v=(s1^2/n1+s2^2/n2)^2/((s1^2/n1)^2/(n1-1)+(s2^2/n2)^2/(n2-1))
    s=sqrt(s1^2/n1+s2^2/n2)
    }
  t0=(barx1-barx2-Delta0)/s

  if(H1=="two")
  {
    t_a=qt(1-alpha/2,df=v)
    cat("H1 is two-tailed. The results are:","\n")
    if(abs(t0)<= t_a){
      cat("\n")
      cat("1. Test statistic t0 is", t0,", t_(v,alpha/2) is", t_a,
          ". Because |t0|<=t_(v,alpha/2), we fail to reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("1. Test statistic t0 is", t0,", t_(v,alpha/2) is", t_a,
          ". Because |t0|>t_(v,alpha/2), we reject H0 at significance level", alpha,"\n")
    }
    Ln=barx1-barx2-t_a*s
    Un=barx1-barx2+t_a*s
    if(Ln<=Delta0 & Delta0<=Un){
      cat("\n")
      cat("2. A", (1-alpha)*100, "% two-tailed confidence interval for the population mean is [", Ln,",",Un,"]","which contains the hypothesized value Delta0=", Delta0, ", so we fail to reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("2. A", (1-alpha)*100, "% two-tailed confidence interval for the population mean is [", Ln,",",Un,"]","which does not contain the hypothesized value Delta0=", Delta0, ", so we reject H0 at significance level", alpha,"\n")
    }
    pv=2*(1-pt(abs(t0),df=v))
    if(pv< alpha){
      cat("\n")
      cat("3. The P-value is", pv ,"which is smaller than alpha=",alpha,", so we reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("3. The P-value is", pv ,"which is not smaller than alpha=",alpha,", so we fail to reject H0 at significance level", alpha,"\n")
    }
  }else if(H1=="left")
  {
    t_a=qt(1-alpha,df=v)
    cat("H1 is left-tailed. The results are:","\n")
    if(t0>= -t_a){
      cat("\n")
      cat("1. Test statistic t0 is", t0,", t_(v,alpha) is", t_a,
          ". Because -t_(v,alpha)<=t0, we fail to reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("1. Test statistic t0 is", t0,", t_(v,alpha) is", t_a,
          ". Because t0<-t_(v,alpha), we reject H0 at significance level", alpha,"\n")
    }
    Un=barx1-barx2+t_a*s
    if(Delta0<=Un){
      cat("\n")
      cat("2. A", (1-alpha)*100, "% one-sided confidence interval for the population mean is (", -Inf,",",Un,"]","which contains the hypothesized value Delta0=", Delta0, ", so we fail to reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("2. A", (1-alpha)*100, "% one-sided confidence interval for the population mean is (", -Inf,",",Un,"]","which does not contain the hypothesized value Delta0=", Delta0, ", so we reject H0 at significance level", alpha,"\n")
    }
    pv=pt(t0,df=v)
    if(pv< alpha){
      cat("\n")
      cat("3. The P-value is", pv ,"which is smaller than alpha=",alpha,", so we reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("3. The P-value is", pv ,"which is not smaller than alpha=",alpha,", so we fail to reject H0 at significance level", alpha,"\n")
    }
  }else if(H1=="right")
  {
    t_a=qt(1-alpha,df=v)
    cat("H1 is right-tailed. The results are:","\n")
    if(t0<= t_a){
      cat("\n")
      cat("1. Test statistic t0 is", t0,", t_(v,alpha) is", t_a,
          ". Because t0<= t_(v,alpha), we fail to reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("1. Test statistic t0 is", t0,", t_(v,alpha) is", t_a,
          ". Because t0>t_(v,alpha), we reject H0 at significance level", alpha,"\n")
    }
    Ln=barx1-barx2-t_a*s
    if(Delta0>=Ln){
      cat("\n")
      cat("2. A", (1-alpha)*100, "% one-sided confidence interval for the population mean is [", Ln,",",Inf,")","which contains the hypothesized value Delta0=", Delta0, ", so we fail to reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("2. A", (1-alpha)*100, "% one-sided confidence interval for the population mean is [", Ln,",",Inf,")","which does not contain the hypothesized value Delta0=", Delta0, ", so we reject H0 at significance level", alpha,"\n")
    }
    pv=1-pt(t0,df=v)
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


Finterval=function(level,sample1,sample2,n1,n2,s1,s2)
{
  if(level>=1|level<=0){return("the confidence level should be between 0 and 1!")}
  if(missing(sample1)==FALSE){
    n1=length(sample1)
    s1=sd(sample1)
  }
  if(missing(sample2)==FALSE){
    n2=length(sample2)
    s2=sd(sample2)
  }
  r=s1^2/s2^2
  alpha=1-level
  fL=qf(alpha/2,n2-1,n1-1)
  fU=qf(1-alpha/2,n2-1,n1-1)
  cat("A", level*100, "% two-sided confidence interval for the ratio between two population variances is [", r*fL ,",",r*fU,"]","\n")
  fL=qf(alpha,n2-1,n1-1)
  fU=qf(1-alpha,n2-1,n1-1)
  cat("A", level*100, "% upper-confidence bound for the ratio between two population variances is", r*fU,"\n")
  cat("A", level*100, "% lower-confidence bound for the ratio between two population variances is", r*fL,"\n")
  cat("\n")
  alpha=1-level
  fL=qf(alpha/2,n2-1,n1-1)
  fU=qf(1-alpha/2,n2-1,n1-1)
  cat("A", level*100, "% two-sided confidence interval for the ratio between two population standard deviations is [", sqrt(r*fL) ,",",sqrt(r*fU),"]","\n")
  fL=qf(alpha,n2-1,n1-1)
  fU=qf(1-alpha,n2-1,n1-1)
  cat("A", level*100, "% upper-confidence bound for the ratio between two population standard deviations is", sqrt(r*fU),"\n")
  cat("A", level*100, "% lower-confidence bound for the ratio between two population standard deviations is", sqrt(r*fL),"\n")
}

Ftest=function(H1="two",alpha,sample1,sample2,n1,s1,n2,s2)
{
  if(alpha>=1|alpha<=0){return("the significance level alpha should be between 0 and 1!")}
  if(missing(sample1)==FALSE){
    n1=length(sample1)
    s1=sd(sample1)
  }
  if(missing(sample2)==FALSE){
    n2=length(sample2)
    s2=sd(sample2)
  }
  F0=s1^2/s2^2

  if(H1=="two")
  {
    fL=qf(alpha/2,n1-1,n2-1)
    fU=qf(1-alpha/2,n1-1,n2-1)

    cat("H1 is two-tailed. The results are:","\n")
    if(F0>=fL & F0<=fU){
      cat("\n")
      cat("1. Test statistic F0 is", F0,", f_(n1-1,n2-1,1-alpha/2) is", fL, ", f_(n1-1,n2-1,alpha/2) is", fU,
          ". Because f_(n1-1,n2-1,1-alpha/2)<=F0<=f_(n1-1,n2-1,alpha/2), we fail to reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("1. Test statistic F0 is", F0,", f_(n1-1,n2-1,1-alpha/2) is", fL, ", f_(n1-1,n2-1,alpha/2) is", fU,
          ". Because f_(n1-1,n2-1,1-alpha/2)<=F0<=f_(n1-1,n2-1,alpha/2) does not hold, we reject H0 at significance level", alpha,"\n")
    }
    fL=qf(alpha/2,n2-1,n1-1)
    fU=qf(1-alpha/2,n2-1,n1-1)
    Ln=F0*fL
    Un=F0*fU
    Delta0=1
    if(Ln<=Delta0 & Delta0<=Un){
      cat("\n")
      cat("2. A", (1-alpha)*100, "% two-tailed confidence interval for the ratio between two population variances is [", Ln,",",Un,"]","which contains 1, so we fail to reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("2. A", (1-alpha)*100, "% two-tailed confidence interval for the ratio between two population variances is [", Ln,",",Un,"]","which does not contain 1, so we reject H0 at significance level", alpha,"\n")
    }
    med=qf(0.5,n1-1,n2-1)
    if(F0<=med){pv=2*pf(F0,n1-1,n2-1)}else{pv=2*(1-pf(F0,n1-1,n2-1))}
    if(pv< alpha){
      cat("\n")
      cat("3. The P-value is", pv ,"which is smaller than alpha=",alpha,", so we reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("3. The P-value is", pv ,"which is not smaller than alpha=",alpha,", so we fail to reject H0 at significance level", alpha,"\n")
    }
  }else if(H1=="left")
  {
    fL=qf(alpha,n1-1,n2-1)

    cat("H1 is left-tailed. The results are:","\n")
    if(F0>=fL){
      cat("\n")
      cat("1. Test statistic F0 is", F0,", f_(n1-1,n2-1,1-alpha) is", fL,
          ". Because f_(n1-1,n2-1,1-alpha)<=F0, we fail to reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("1. Test statistic F0 is", F0,", f_(n1-1,n2-1,1-alpha) is", fL,
          ". Because F0<=f_(n1-1,n2-1,1-alpha), we reject H0 at significance level", alpha,"\n")
    }
    fU=qf(1-alpha,n2-1,n1-1)
    Un=F0*fU
    Delta0=1
    if(Delta0<=Un){
      cat("\n")
      cat("2. A", (1-alpha)*100, "% two-tailed confidence interval for the ratio between two population variances is [", -Inf,",",Un,"]","which contains 1, so we fail to reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("2. A", (1-alpha)*100, "% two-tailed confidence interval for the ratio between two population variances is [", -Inf,",",Un,"]","which does not contain 1, so we reject H0 at significance level", alpha,"\n")
    }
    pv=pf(F0,n1-1,n2-1)
    if(pv< alpha){
      cat("\n")
      cat("3. The P-value is", pv ,"which is smaller than alpha=",alpha,", so we reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("3. The P-value is", pv ,"which is not smaller than alpha=",alpha,", so we fail to reject H0 at significance level", alpha,"\n")
    }
  }else if(H1=="right")
  {
    fU=qf(1-alpha,n1-1,n2-1)

    cat("H1 is right-tailed. The results are:","\n")
    if(F0<=fU){
      cat("\n")
      cat("1. Test statistic F0 is", F0,", f_(n1-1,n2-1,alpha) is", fU,
          ". Because F0<=f_(n1-1,n2-1,alpha), we fail to reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("1. Test statistic F0 is", F0,", f_(n1-1,n2-1,alpha) is", fU,
          ". Because F0>f_(n1-1,n2-1,alpha), we reject H0 at significance level", alpha,"\n")
    }
    fL=qf(alpha,n2-1,n1-1)
    Ln=F0*fL
    Delta0=1
    if(Ln<=Delta0){
      cat("\n")
      cat("2. A", (1-alpha)*100, "% two-tailed confidence interval for the ratio between two population variances is [", Ln,",",Inf,"]","which contains 1, so we fail to reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("2. A", (1-alpha)*100, "% two-tailed confidence interval for the ratio between two population variances is [", Ln,",",Inf,"]","which does not contain 1, so we reject H0 at significance level", alpha,"\n")
    }
    pv=1-pf(F0,n1-1,n2-1)
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



twosample.Propinterval=function(level, n1,n2,X1,X2)
{
  if(level>=1|level<=0){return("the confidence level should be between 0 and 1!")}
  if(missing(n1)==TRUE | missing(X1)==TRUE|missing(n2)==TRUE | missing(X2)==TRUE)
  {
    return("No samples were provided, please input the sample sizes/the values of X1 and X2")
  }
  if(X1<5 | n1-X1<5|X2<5|n2-X2<5){cat("Warning: Because at least one of X1,n1-X1,X2, and n2-X2, is less than 5. The CLT might not work well, and the following results should be used with caution ")}
  p1=X1/n1
  p2=X2/n2
  alpha=1-level
  z_a=qnorm(1-alpha/2)
  s=sqrt(p1*(1-p1)/n1+p2*(1-p2)/n2)
  cat("A large-sample confidence interval for the difference in population proportions with confidence level of approximately", level*100, "% is [", p1-p2-z_a*s,",",p1-p2+z_a*s,"]","\n")
  z_a=qnorm(1-alpha)
  cat("A large-sample upper-confidence bound for the difference in population proportions with confidence level of approximately", level*100, "% is", p1-p2+z_a*s,"\n")
  cat("A large-sample lower-confidence bound for the difference in population proportions with confidence level of approximately", level*100, "% is", p1-p2-z_a*s,"\n")
}

twosample.Proptest=function(H1="two",alpha,n1,n2,X1,X2)
{
  if(alpha>=1|alpha<=0){return("the significance level alpha should be between 0 and 1!")}
  if(missing(n1)==TRUE | missing(X1)==TRUE|missing(n2)==TRUE | missing(X2)==TRUE)
  {
    return("No samples were provided, please input the sample sizes/the values of X1 and X2")
  }
  if(X1<5 | n1-X1<5|X2<5|n2-X2<5){cat("Warning: Because at least one of X1,n1-X1,X2, and n2-X2, is less than 5. The CLT might not work well, and the following results should be used with caution ")}
  p1=X1/n1
  p2=X2/n2
  p=(X1+X2)/(n1+n2)

  z0=(p1-p2)/sqrt(p*(1-p)*(1/n1+1/n2))
  s=sqrt(p1*(1-p1)/n1+p2*(1-p2)/n2)

  if(H1=="two"){
    z_a=qnorm(1-alpha/2)
    cat("H1 is two-tailed. The results are:","\n")
    if(abs(z0)<=z_a){
      cat("\n")
      cat("1. Test statistic z0 is", z0,", z_(alpha/2) is", z_a,
          ". Because |z0|<=z_(alpha/2), we fail to reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("1. Test statistic z0 is", z0,", z_(alpha/2) is", z_a,
          ". Because |z0|>z_(alpha/2), so we reject H0 at significance level", alpha,"\n")
    }
    Ln=p1-p2-z_a*s
    Un=p1-p2+z_a*s
    if(Ln<=0 & 0<=Un){
      cat("\n")
      cat("2. A", (1-alpha)*100, "% two-tailed confidence interval for the difference in population proportions is [", Ln,",",Un,"]","which contains 0, so we fail to reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("2. A", (1-alpha)*100, "% two-tailed confidence interval for the difference in population proportions is [", Ln,",",Un,"]","which does not contain 0, so we reject H0 at significance level", alpha,"\n")
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
    cat("H1 is left-tailed. The results are:","\n")
    if(z0>= -z_a){
      cat("\n")
      cat("1. Test statistic z0 is", z0,", z_alpha is", z_a,
          ". Because -z_alpha<=z0, we fail to reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("1. Test statistic z0 is", z0,", z_alpha is", z_a,
          ". Because z0<-z_alpha, so we reject H0 at significance level", alpha,"\n")
    }
    Un=p1-p2+z_a*s
    if(0<=Un){
      cat("\n")
      cat("2. A", (1-alpha)*100, "% one-sided confidence interval for the difference in population proportions is (", -Inf,",",Un,"]","which contains 0, so we fail to reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("2. A", (1-alpha)*100, "% one-sided confidence interval for the difference in population proportions is (", -Inf,",",Un,"]","which does not contain 0, so we reject H0 at significance level", alpha,"\n")
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
    cat("H1 is right-tailed. The results are:","\n")
    if(z0<= z_a){
      cat("\n")
      cat("1. Test statistic z0 is", z0,", z_alpha is", z_a,
          ". Because z0<= z_alpha, we fail to reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("1. Test statistic z0 is", z0,", z_alpha is", z_a,
          ". Because z0>z_alpha, so we reject H0 at significance level", alpha,"\n")
    }
    Ln=p1-p2+z_a*s
    if(Delta0>=Ln){
      cat("\n")
      cat("2. A", (1-alpha)*100, "% one-sided confidence interval for the difference in population proportions is [", Ln,",",Inf,")","which contains 0, so we fail to reject H0 at significance level", alpha,"\n")
    }else{
      cat("\n")
      cat("2. A", (1-alpha)*100, "% one-sided confidence interval for the difference in population proportions is [", Ln,",",Inf,")","which does not contain 0, so we reject H0 at significance level", alpha,"\n")
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
