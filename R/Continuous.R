f=function(x)
{
  res=I(x>12.5)*20*exp(-20*(x-12.5))
  return(res)
}

continuous.plotpdf=function(f,a,b)
{
  x=seq(a,b,length=100)
  fx=sapply(x,f)
  par(mar=c(4,4,.5,.5))
  plot(x,fx,ylim=c(0,max(fx)*1.1),xlab="x",ylab="f(x)",type="l")
}

continuous.plotcdf=function(f,a,b)
{
  cdf=function(a){return(integrate(f,-Inf,a))}
  x=seq(a,b,length=100)
  Fx=sapply(x,cdf)
  par(mar=c(4,4,.5,.5))
  plot(x,fx,ylim=c(0,max(fx)*1.1),xlab="x",ylab="f(x)",type="l")
}
