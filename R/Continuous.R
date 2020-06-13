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
