###################################
#Continuous Uniform Distribution
uniform.summary=function(a,b,plotpdf=TRUE,plotcdf=TRUE)
{
  if(a>=b){return("a must be smaller than b")}
  if(a==-Inf |b==Inf|a==Inf|b==-Inf){return("a and b must be finite")}
  mu=(a+b)/2
  sigma2=(b-a)^2/12
  sigma=sqrt(sigma2)

  if(plotpdf==TRUE & plotcdf==FALSE)
  {
    s=seq(a,b,length=100)
    plot(s,dunif(s,a,b),xlab="x",ylab="f(x)",type="l",ylim=c(0,1.1/(b-a)),xlim=c(a-(b-a)/10,b+(b-a)/10))
    lines(seq(a-(b-a)/10,a,length=100),rep(0,100))
    lines(seq(b,b+(b-a)/10,length=100),rep(0,100))
    segments(a,0,a,1/(b-a),lty=3)
    segments(b,0,b,1/(b-a),lty=3)
  }
  if(plotpdf==FALSE & plotcdf==TRUE)
  {
    s=seq(a-(b-a)/10,b+(b-a)/10,length=100)
    plot(s,punif(s,a,b),xlab="x",ylab="F(x)",type="l",ylim=c(0,1))
  }
  if(plotpdf==TRUE & plotcdf==TRUE)
  {
    par(mfrow=c(2,1))
    par(mar=c(4,4,.5,.1))
    s=seq(a,b,length=100)
    plot(s,dunif(s,a,b),xlab="x",ylab="f(x)",type="l",ylim=c(0,1.1/(b-a)),xlim=c(a-(b-a)/10,b+(b-a)/10))
    lines(seq(a-(b-a)/10,a,length=100),rep(0,100))
    lines(seq(b,b+(b-a)/10,length=100),rep(0,100))
    segments(a,0,a,1/(b-a),lty=3)
    segments(b,0,b,1/(b-a),lty=3)
    s=seq(a-(b-a)/10,b+(b-a)/10,length=100)
    plot(s,punif(s,a,b),xlab="x",ylab="F(x)",type="l",ylim=c(0,1))
  }
  par(mfrow=c(1,1))
  return(list(mean=mu,variance=sigma2,standard.deviation=sigma))
}

uniform.prob=function(a,b,lb,ub)
{
  if(a>=b){return("a must be smaller than b")}
  if(a==-Inf |b==Inf|a==Inf|b==-Inf){return("a and b must be finite")}
  if(ub<lb){return("lb must be smaller than ub!")}
  return(punif(ub,a,b)-punif(lb,a,b))
}

uniform.quantile=function(a,b,q)
{
  if(a>=b){return("a must be smaller than b")}
  if(a==-Inf |b==Inf|a==Inf|b==-Inf){return("a and b must be finite")}
  if(q<=0|q>=1){return("q must be between 0 and 1")}
  return(qunif(q,a,b))
}

###################################
#Normal Distribution
normal.summary=function(mu,sigma,plotpdf=TRUE,plotcdf=TRUE)
{
  if(abs(mu)==Inf | abs(sigma)==Inf){return("mu and sigma must be finite")}
  if(sigma<=0){return("sigma must be positive")}
  mu=mu
  sigma2=sigma^2
  sigma=sigma

  if(plotpdf==TRUE & plotcdf==FALSE)
  {
    s=seq(mu-4*sigma,mu+4*sigma,length=200)
    plot(s,dnorm(s,mu,sigma),xlab="x",ylab="f(x)",type="l")
  }
  if(plotpdf==FALSE & plotcdf==TRUE)
  {
    s=seq(mu-4*sigma,mu+4*sigma,length=200)
    plot(s,pnorm(s,mu,sigma),xlab="x",ylab="F(x)",type="l",ylim=c(0,1))
  }
  if(plotpdf==TRUE & plotcdf==TRUE)
  {
    par(mfrow=c(2,1))
    par(mar=c(4,4,.5,.1))
    s=seq(mu-4*sigma,mu+4*sigma,length=200)
    plot(s,dnorm(s,mu,sigma),xlab="x",ylab="f(x)",type="l")
    plot(s,pnorm(s,mu,sigma),xlab="x",ylab="F(x)",type="l",ylim=c(0,1))
  }
  par(mfrow=c(1,1))
  return(list(mean=mu,variance=sigma2,standard.deviation=sigma))
}

normal.prob=function(mu,sigma,lb,ub)
{
  if(abs(mu)==Inf | abs(sigma)==Inf){return("mu and sigma must be finite")}
  if(sigma<=0){return("sigma must be positive")}
  if(ub<lb){return("lb must be smaller than ub!")}
  return(pnorm(ub,mu,sigma)-pnorm(lb,mu,sigma))
}

normal.quantile=function(mu,sigma,q)
{
  if(abs(mu)==Inf | abs(sigma)==Inf){return("mu and sigma must be finite")}
  if(sigma<=0){return("sigma must be positive")}
  if(q<=0|q>=1){return("q must be between 0 and 1")}
  return(qnorm(q,mu,sigma))
}


###################################
#Exponential Distribution
exponential.summary=function(lambda,plotpdf=TRUE,plotcdf=TRUE)
{
  if(abs(lambda)==Inf | lambda<=0){return("lambda must be a finite positive number")}
  mu=1/lambda
  sigma2=1/lambda^2
  sigma=1/lambda

  if(plotpdf==TRUE & plotcdf==FALSE)
  {
    s=seq(0,qexp(0.999,rate=lambda),length=100)
    plot(s,lambda*exp(-lambda*s),xlab="x",ylab="f(x)",type="l")
  }
  if(plotpdf==FALSE & plotcdf==TRUE)
  {
    s=seq(0,qexp(0.999,rate=lambda),length=100)
    plot(s,1-exp(-lambda*s),xlab="x",ylab="F(x)",type="l",ylim=c(0,1))
  }
  if(plotpdf==TRUE & plotcdf==TRUE)
  {
    par(mfrow=c(2,1))
    par(mar=c(4,4,.5,.1))
    s=seq(0,qexp(0.999,rate=lambda),length=100)
    plot(s,lambda*exp(-lambda*s),xlab="x",ylab="f(x)",type="l")
    plot(s,1-exp(-lambda*s),xlab="x",ylab="F(x)",type="l",ylim=c(0,1))
  }
  par(mfrow=c(1,1))
  return(list(mean=mu,variance=sigma2,standard.deviation=sigma))
}

exponential.prob=function(lambda,lb,ub)
{
  if(abs(lambda)==Inf | lambda<=0){return("lambda must be a finite positive number")}
  if(ub<lb){return("lb must be smaller than ub!")}
  if(lb>0){return(exp(-lambda*lb)-exp(-lambda*ub))}
  if(ub<=0){return(0)}
  if(lb<=0 & ub>0){return(1-exp(-lambda*ub))}
}

exponential.quantile=function(lambda,q)
{
  if(abs(lambda)==Inf | lambda<=0){return("lambda must be a finite positive number")}
  if(q<=0|q>=1){return("q must be between 0 and 1")}
  return(-log(1-q)/lambda)
}

###################################
#Exponential Distribution
gamma.summary=function(r,lambda,plotpdf=TRUE,plotcdf=TRUE)
{
  if(abs(lambda)==Inf | lambda<=0){return("lambda must be a finite positive number")}
  if(abs(r)==Inf | r<=0){return("r must be a finite positive number")}

  mu=r/lambda
  sigma2=r/lambda^2
  sigma=sqrt(sigma2)

  if(plotpdf==TRUE & plotcdf==FALSE)
  {
    s=seq(0,qgamma(0.999,shape=r,scale=lambda),length=100)
    plot(s,dgamma(s,shape=r,scale=1/lambda),xlab="x",ylab="f(x)",type="l")
  }
  if(plotpdf==FALSE & plotcdf==TRUE)
  {
    s=seq(0,qgamma(0.999,shape=r,scale=lambda),length=100)
    plot(s,pgamma(s,shape=r,scale=1/lambda),xlab="x",ylab="F(x)",type="l",ylim=c(0,1))
  }
  if(plotpdf==TRUE & plotcdf==TRUE)
  {
    par(mfrow=c(2,1))
    par(mar=c(4,4,.5,.1))
    s=seq(0,qgamma(0.999,shape=r,scale=1/lambda),length=100)
    plot(s,dgamma(s,shape=r,scale=1/lambda),xlab="x",ylab="f(x)",type="l")
    plot(s,pgamma(s,shape=r,scale=1/lambda),xlab="x",ylab="F(x)",type="l",ylim=c(0,1))
  }
  par(mfrow=c(1,1))
  return(list(mean=mu,variance=sigma2,standard.deviation=sigma))
}

gamma.prob=function(r,lambda,lb,ub)
{
  if(abs(lambda)==Inf | lambda<=0){return("lambda must be a finite positive number")}
  if(abs(r)==Inf | r<=0){return("r must be a finite positive number")}
  if(ub<lb){return("lb must be smaller than ub!")}
  return(pgamma(ub,shape=r,scale=1/lambda)-pgamma(lb,shape=r,scale=1/lambda))
}

gamma.quantile=function(r,lambda,q)
{
  if(abs(lambda)==Inf | lambda<=0){return("lambda must be a finite positive number")}
  if(abs(r)==Inf | r<=0){return("r must be a finite positive number")}
  if(q<=0|q>=1){return("q must be between 0 and 1")}
  return(qgamma(q,shape=r,scale=1/lambda))
}


###################################
#Weibull Distribution
weibull.summary=function(beta,delta,plotpdf=TRUE,plotcdf=TRUE)
{
  if(abs(delta)==Inf | delta<=0){return("delta must be a finite positive number")}
  if(abs(beta)==Inf | beta<=0){return("beta must be a finite positive number")}

  mu=delta*gamma(1+1/beta)
  sigma2=delta^2*(gamma(1+2/beta)-(gamma(1+1/beta))^2)
  sigma=sqrt(sigma2)

  if(plotpdf==TRUE & plotcdf==FALSE)
  {
    s=seq(0,qweibull(0.999,shape=beta,scale=delta),length=100)
    plot(s,dweibull(s,shape=beta,scale=delta),xlab="x",ylab="f(x)",type="l")
  }
  if(plotpdf==FALSE & plotcdf==TRUE)
  {
    s=seq(0,qweibull(0.999,shape=beta,scale=delta),length=100)
    plot(s,pweibull(s,shape=beta,scale=delta),xlab="x",ylab="F(x)",type="l",ylim=c(0,1))
  }
  if(plotpdf==TRUE & plotcdf==TRUE)
  {
    par(mfrow=c(2,1))
    par(mar=c(4,4,.5,.1))
    s=seq(0,qweibull(0.999,shape=beta,scale=delta),length=100)
    plot(s,dweibull(s,shape=beta,scale=delta),xlab="x",ylab="f(x)",type="l")
    plot(s,pweibull(s,shape=beta,scale=delta),xlab="x",ylab="F(x)",type="l",ylim=c(0,1))
  }
  par(mfrow=c(1,1))
  return(list(mean=mu,variance=sigma2,standard.deviation=sigma))
}

weibull.prob=function(beta,delta,lb,ub)
{
  if(abs(delta)==Inf | delta<=0){return("delta must be a finite positive number")}
  if(abs(beta)==Inf | beta<=0){return("beta must be a finite positive number")}
  if(ub<lb){return("lb must be smaller than ub!")}
  return(pweibull(ub,shape=beta,scale=delta)-pweibull(lb,shape=beta,scale=delta))
}

weibull.quantile=function(beta,delta,q)
{
  if(abs(delta)==Inf | delta<=0){return("delta must be a finite positive number")}
  if(abs(beta)==Inf | beta<=0){return("beta must be a finite positive number")}
  if(q<=0|q>=1){return("q must be between 0 and 1")}
  return(qweibull(q,shape=beta,scale=delta))
}


###################################
#Lognormal Distribution
lognormal.summary=function(theta,omega,plotpdf=TRUE,plotcdf=TRUE)
{
  if(abs(omega)==Inf | omega<=0){return("omega must be a finite positive number")}
  if(abs(theta)==Inf){return("theta must be a finite number")}

  mu=exp(theta+omega^2/2)
  sigma2=exp(2*theta+omega^2)*(exp(omega^2)-1)
  sigma=sqrt(sigma2)

  if(plotpdf==TRUE & plotcdf==FALSE)
  {
    s=seq(0,qlnorm(0.99,theta,omega),length=1000)
    plot(s,dlnorm(s,theta,omega),xlab="x",ylab="f(x)",type="l")
  }
  if(plotpdf==FALSE & plotcdf==TRUE)
  {
    s=seq(0,qlnorm(0.99,theta,omega),length=1000)
    plot(s,plnorm(s,theta,omega),xlab="x",ylab="F(x)",type="l",ylim=c(0,1))
  }
  if(plotpdf==TRUE & plotcdf==TRUE)
  {
    par(mfrow=c(2,1))
    par(mar=c(4,4,.5,.1))
    s=seq(0,qlnorm(0.99,theta,omega),length=1000)
    plot(s,dlnorm(s,theta,omega),xlab="x",ylab="f(x)",type="l")
    plot(s,plnorm(s,theta,omega),xlab="x",ylab="F(x)",type="l",ylim=c(0,1))
  }
  par(mfrow=c(1,1))
  return(list(mean=mu,variance=sigma2,standard.deviation=sigma))
}

lognormal.prob=function(theta,omega,lb,ub)
{
  if(abs(omega)==Inf | omega<=0){return("omega must be a finite positive number")}
  if(abs(theta)==Inf){return("theta must be a finite number")}
  if(ub<lb){return("lb must be smaller than ub!")}
  return(plnorm(ub,theta,omega)-plnorm(lb,theta,omega))
}

lognormal.quantile=function(theta,omega,q)
{
  if(abs(omega)==Inf | omega<=0){return("omega must be a finite positive number")}
  if(abs(theta)==Inf){return("theta must be a finite number")}
  if(q<=0|q>=1){return("q must be between 0 and 1")}
  return(qlnorm(q,theta,omega))
}


###################################
#Beta Distribution
beta.summary=function(alpha,beta,plotpdf=TRUE,plotcdf=TRUE)
{
  if(abs(alpha)==Inf | alpha<=0){return("alpha must be a finite positive number")}
  if(abs(beta)==Inf | beta<=0){return("beta must be a finite positive number")}

  mu=alpha/(alpha+beta)
  sigma2=alpha*beta/(alpha+beta)^2/(alpha+beta+1)
  sigma=sqrt(sigma2)

  if(plotpdf==TRUE & plotcdf==FALSE)
  {
    s=seq(0,1,length=1000)
    y=dbeta(s,alpha,beta)
    plot(s,y,xlab="x",ylab="f(x)",type="l",xlim=c(-.15,1.15),ylim=c(0,max(y)+0.02))
    lines(seq(-.15,0,length=100),rep(0,100))
    lines(seq(1,1.15,length=100),rep(0,100))
    segments(0,0,0,y[1],lty=3)
    segments(1,0,1,y[length(y)],lty=3)
  }
  if(plotpdf==FALSE & plotcdf==TRUE)
  {
    s=seq(-.15,1.15,length=1500)
    plot(s,pbeta(s,alpha,beta),xlab="x",ylab="F(x)",type="l",ylim=c(0,1))
  }
  if(plotpdf==TRUE & plotcdf==TRUE)
  {
    par(mfrow=c(2,1))
    par(mar=c(4,4,.5,.1))
    s=seq(0,1,length=1000)
    y=dbeta(s,alpha,beta)
    plot(s,y,xlab="x",ylab="f(x)",type="l",xlim=c(-.15,1.15),ylim=c(0,max(y)+0.02))
    lines(seq(-.15,0,length=100),rep(0,100))
    lines(seq(1,1.15,length=100),rep(0,100))
    segments(0,0,0,y[1],lty=3)
    segments(1,0,1,y[length(y)],lty=3)
    s=seq(-.15,1.15,length=1500)
    plot(s,pbeta(s,alpha,beta),xlab="x",ylab="F(x)",type="l",ylim=c(0,1))
  }
  par(mfrow=c(1,1))
  return(list(mean=mu,variance=sigma2,standard.deviation=sigma))
}

beta.prob=function(alpha,beta,lb,ub)
{
  if(abs(alpha)==Inf | alpha<=0){return("alpha must be a finite positive number")}
  if(abs(beta)==Inf | beta<=0){return("beta must be a finite positive number")}
  if(ub<lb){return("lb must be smaller than ub!")}
  return(pbeta(ub,alpha,beta)-pbeta(lb,alpha,beta))
}

beta.quantile=function(alpha,beta,q)
{
  if(abs(alpha)==Inf | alpha<=0){return("alpha must be a finite positive number")}
  if(abs(beta)==Inf | beta<=0){return("beta must be a finite positive number")}
  if(q<=0|q>=1){return("q must be between 0 and 1")}
  return(qbeta(q,alpha,beta))
}
