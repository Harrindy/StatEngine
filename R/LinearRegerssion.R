lm.est=function(fit)
{
  p=length(fit$coefficients)
  name=c()
  for(k in 1:p)
  {
    name=c(name,paste("beta",k-1,sep=""))
  }
  name=c(name,"sigma")
  Estimate=c(fit$coefficients,sqrt(sum(fit$residuals^2)/fit$df.residual))
  return(data.frame(Estimate,row.names = name))
}

lm.coef.interval=function(fit,level=0.95)
{
  p=length(fit$coefficients)
  name=c()
  for(k in 1:p)
  {
    name=c(name,paste("beta",k-1,sep=""))
  }
  se=sqrt(diag(vcov(fit)))
  DF=fit$df.residual
  alpha=1-level
  tt=qt(1-alpha/2,DF)
  LS.est=fit$coefficients
  CI.lb=LS.est-tt*se
  CI.ub=LS.est+tt*se
  cat("Two-sided", level*100, "% confidence intervals of regression coefficients are \n")
  print(data.frame(cbind(CI.lb,CI.ub),row.names = name))
  cat("\n")
  tt=qt(1-alpha,DF)
  upper.bound=LS.est+tt*se
  lower.bound=LS.est-tt*se
  cat("One-sided",level*100, "% (lower and upper) confidence bounds are \n")
  print(data.frame(cbind(lower.bound,upper.bound),row.names = name))
}

lm.coef.test=function(fit,alpha=0.05,H1="two",hypo.beta=rep(0,length(fit$coefficients)),round=6)
{
  p=length(fit$coefficients)
  name=c()
  for(k in 1:p)
  {
    name=c(name,paste("beta",k-1,sep=""))
  }
  if(length(hypo.beta)!=p){return("the length of beta0 should be p")}

  LS.est=fit$coefficients
  se=sqrt(diag(vcov(fit)))
  T0=(LS.est-hypo.beta)/se
  if(H1=="two"){
    T0.abs=abs(T0)
    tt=qt(1-alpha/2,df=fit$df.residual)
    t_alpha.2=tt
    P_value=2*(1-pt(T0.abs,df=fit$df.residual))
    esti.beta=fit$coefficients
    CI.lb=esti.beta-tt*se
    CI.ub=esti.beta+tt*se
    res=P_value<alpha
    Reject=rep("No",p)
    for(k in 1:p)
    {
      if(res[k]==TRUE){Reject[k]="Yes"}
    }
    res=round(cbind(P_value,t_alpha.2,T0.abs,CI.lb,hypo.beta,CI.ub),digits=round)
    print(data.frame(cbind(res,Reject),row.names=name))
  }else if(H1=="left"){
    tt=-qt(1-alpha,df=fit$df.residual)
    t_alpha=tt
    P_value=pt(T0,df=fit$df.residual)
    esti.beta=fit$coefficients
    CI.lb=-Inf
    CI.ub=esti.beta+tt*se
    res=P_value<alpha
    Reject=rep("No",p)
    for(k in 1:p)
    {
      if(res[k]==TRUE){Reject[k]="Yes"}
    }
    res=round(cbind(P_value,T0,t_alpha,CI.lb,hypo.beta,CI.ub),digits=round)
    print(data.frame(cbind(res,Reject),row.names=name))
  }else if(H1=="right"){
    tt=qt(1-alpha,df=fit$df.residual)
    t_alpha=tt
    P_value=1-pt(T0,df=fit$df.residual)
    esti.beta=fit$coefficients
    CI.ub=Inf
    CI.lb=esti.beta-tt*se
    res=P_value<alpha
    Reject=rep("No",p)
    for(k in 1:p)
    {
      if(res[k]==TRUE){Reject[k]="Yes"}
    }
    res=round(cbind(P_value,t_alpha,T0,CI.lb,hypo.beta,CI.ub),digits=round)
    print(data.frame(cbind(res,Reject),row.names=name))
  }
}

lm.partialFtest=function(fit.H0,fit.ALL,alpha=0.05,round=6)
{
  anova.h0=anova(fit.H0)
  anova.h1=anova(fit.ALL)
  r=length(fit.ALL$coefficients)-length(fit.H0$coefficients)
  DF=fit.ALL$df.residual
  F0=(tail(anova.h0[[2]],1)-tail(anova.h1[[2]],1))/r/tail(anova.h1[[3]],1)
  f_alpha=qf(1-alpha,r,DF)
  P_value=1-pf(F0,r,DF)
  Reject="No"
  if(F0>f_alpha){Reject="Yes"}
  res=round(cbind(F0,f_alpha,P_value),digits=round)
  print(data.frame(cbind(res,Reject),row.names = "PartialFtest"))
}

lm.modelcheck=function(fit)
{
  require(car)
  par(mfrow=c(1,3))
  par(mar=c(4,4,2,.5))
  e=fit$residuals
  qqPlot(e,xlab="normal quantile",ylab="Observations",main="QQ-plot of ordinary residuals")
  d=e/(sum(e^2)/fit$df.residual)
  plot(d,ylim=c(min(-4,min(d)),max(4,max(d))),main="Standarized residuals")
  lines(0:(length(d)+1),rep(3,length(d)+2),col="blue",lty=5,lwd=2)
  lines(0:(length(d)+1),rep(-3,length(d)+2),col="blue",lty=5,lwd=2)
  D=cooks.distance(fit)
  plot(D,ylim=c(0,max(1,max(D))),main="Cook's Distances")
  lines(0:(length(D)+1),rep(1,length(D)+2),col="blue",lty=5,lwd=2)
  par(mfrow=c(1,1))
  cat("VIFs are: \n")
  v=vif(fit)
  print(v)
  #return(list(ordinary.residuals=e,standarized.residuals=d,Cook.distance=D,VIF=v))
}
