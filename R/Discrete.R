
discrete.plotpdf=function(x,fx)
{
  #if(abs(sum(fx)-1)>1e-10){return("The summation of entries in fx must be 1!")}
  if(min(fx)<0){return("Entries in fx must be non-negative!")}
  if(length(unique(x))<length(x)){return("Entries in x must be unique!")}
  x=x[fx>0]
  fx=fx[fx>0]
  fx=fx[order(x)]
  x=x[order(x)]
  k=length(x)
  par(mar=c(4,4,.5,.5))
  plot(x,fx,xlim=c(min(x)-(x[2]-x[1])/2,max(x)+(x[length(x)]-x[length(x)-1])/2),ylim=c(0,max(fx)*1.2),xlab="x",ylab="f(x)")
  lines(seq(min(x)-(x[2]-x[1])/2,max(x)+(x[length(x)]-x[length(x)-1])/2,length=100),rep(0,100))
  points(x,fx,pch=16)
  for(i in 1:k)
  {
    segments(x[i],0,x[i],fx[i],lwd=2)
  }
}

discrete.plotcdf=function(x,fx)
{
  #if(abs(sum(fx)-1)>1e-10){return("The summation of entries in fx must be 1!")}
  if(min(fx)<0){return("Entries in fx must be non-negative!")}
  if(length(unique(x))<length(x)){return("Entries in x must be unique!")}
  x=x[fx>0]
  fx=fx[fx>0]
  fx=fx[order(x)]
  x=x[order(x)]
  k=length(x)
  Fx=c(cumsum(fx))

  par(mar=c(4,4,.5,.5))
  plot(x,Fx,xlim=c(min(x)-(x[2]-x[1])/2,max(x)+(x[length(x)]-x[length(x)-1])/2),ylim=c(0,1),xlab="x",ylab="F(x)")
  lines(seq(min(x)-(x[2]-x[1])/2,min(x),length=100),rep(0,100),lwd=2)
  points(x[1],0,pch=1)
  if(abs(sum(fx)-1)<=1e-10){lines(seq(max(x),max(x)+(x[length(x)]-x[length(x)-1])/2,length=100),rep(1,100),lwd=2)}
  points(x,Fx,pch=16)
  for(i in 1:(k-1))
  {
    lines(seq(x[i],x[i+1],length=100),rep(Fx[i],100),lwd=2)
    points(x[i+1],Fx[i],pch=1)
  }
}

discrete.summary=function(x,fx,plotpdf=TRUE,plotcdf=TRUE)
{
  if(abs(sum(fx)-1)>1e-10){return("The summation of entries in fx must be 1!")}
  if(min(fx)<0){return("Entries in fx must be non-negative!")}
  if(length(unique(x))<length(x)){return("Entries in x must be unique!")}
  x=x[fx>0]
  fx=fx[fx>0]
  mu=sum(x*fx)
  sigma2=sum(x^2*fx)-mu^2
  sigma=sqrt(sigma2)

  if(plotpdf==TRUE & plotcdf==FALSE)
  {
    discrete.plotpdf(x,fx)
  }
  if(plotpdf==FALSE & plotcdf==TRUE)
  {
    discrete.plotcdf(x,fx)
  }
  if(plotpdf==TRUE & plotcdf==TRUE)
  {
    par(mfrow=c(2,1))
    discrete.plotpdf(x,fx)
    discrete.plotcdf(x,fx)
  }
  par(mfrow=c(1,1))
  return(list(mean=mu,variance=sigma2,standard.deviation=sigma))
}

discrete.prob=function(x,fx,lb,ub=lb,inclusive="both")
{
  if(ub<lb){return("ub cannot be smaller than lb!")}
  if(abs(sum(fx)-1)>1e-10){return("The summation of entries in fx must be 1!")}
  if(min(fx)<0){return("Entries in fx must be non-negative!")}
  if(length(unique(x))<length(x)){return("Entries in x must be unique!")}
  x=x[fx>0]
  fx=fx[fx>0]

  if(inclusive=="both")
  {
    return(sum(fx[lb<=x &x<=ub]))
  }else if(inclusive=="left"){
    return(sum(fx[lb<=x &x<ub]))
  }else if(inclusive=="right"){
    return(sum(fx[lb<x &x<=ub]))
  }else if(inclusive=="none"){
    return(sum(fx[lb<x &x<ub]))
  }else{
    return("inclusive must be one of these: both, left, right, none")
  }
}


###################################
#Discrete Uniform Distribution
duniform.summary=function(range,plotpdf=TRUE,plotcdf=TRUE)
{
  x=range
  if(length(unique(x))<length(x)){return("Entries in x must be unique!")}
  fx=rep(1/length(x),length(x))
  mu=sum(x*fx)
  sigma2=sum(x^2*fx)-mu^2
  sigma=sqrt(sigma2)

  if(plotpdf==TRUE & plotcdf==FALSE)
  {
    discrete.plotpdf(x,fx)
  }
  if(plotpdf==FALSE & plotcdf==TRUE)
  {
    discrete.plotcdf(x,fx)
  }
  if(plotpdf==TRUE & plotcdf==TRUE)
  {
    par(mfrow=c(2,1))
    discrete.plotpdf(x,fx)
    discrete.plotcdf(x,fx)
  }
  par(mfrow=c(1,1))
  return(list(mean=mu,variance=sigma2,standard.deviation=sigma))
}

duniform.prob=function(range,lb,ub=lb,inclusive="both")
{
  if(ub<lb){return("ub cannot be smaller than lb!")}
  x=range
  if(length(unique(x))<length(x)){return("Entries in x must be unique!")}
  fx=rep(1/length(x),length(x))
  return(discrete.prob(x,fx,lb,ub,inclusive=inclusive))
}



###################################
#Binomial distribution
binomial.summary=function(n,p,plotpdf=TRUE,plotcdf=TRUE)
{
  if((as.integer(n)-n)!=0 | n<=0){return("n must be a positive integer")}
  if(p<0 | p>1){return("p must be 0<p<1")}
  if(p==0){
    x=0;fx=1
    mu=0
    sigma2=0
    sigma=sqrt(sigma2)
  }else if(p==1){
    x=n;fx=1
    mu=n
    sigma2=0
    sigma=sqrt(sigma2)
  }else{
  x=0:n
  fx=dbinom(x,n,p)
  mu=n*p
  sigma2=n*p*(1-p)
  sigma=sqrt(sigma2)
  }

  if(plotpdf==TRUE & plotcdf==FALSE)
  {
    discrete.plotpdf(x,fx)
  }
  if(plotpdf==FALSE & plotcdf==TRUE)
  {
    discrete.plotcdf(x,fx)
  }
  if(plotpdf==TRUE & plotcdf==TRUE)
  {
    par(mfrow=c(2,1))
    discrete.plotpdf(x,fx)
    discrete.plotcdf(x,fx)
  }
  par(mfrow=c(1,1))
  return(list(mean=mu,variance=sigma2,standard.deviation=sigma))
}

binomial.prob=function(n,p,lb,ub=lb,inclusive="both")
{
  if(ub<lb){return("ub cannot be smaller than lb!")}
  if((as.integer(n)-n)!=0 | n<=0){return("n must be a positive integer")}
  if(p<0 | p>1){return("p must be 0<p<1")}
  if(p==0){
    x=0;fx=1
  }else if(p==1){
    x=n;fx=1
  }else{
    x=0:n
    fx=dbinom(x,n,p)
  }
  return(discrete.prob(x,fx,lb,ub,inclusive=inclusive))
}


###################################
#Negative Binomial distribution
negbinom.summary=function(r,p,plotpdf=TRUE,plotcdf=TRUE)
{
  if((as.integer(r)-r)!=0 | r<=0){return("r must be a positive integer")}
  if(p<=0 | p>1){return("p must be 0<p<1")}

  mu=r/p
  sigma2=r*(1-p)/(p^2)
  sigma=sqrt(sigma2)

  x=r:max(r+10,ceiling(mu)+10)
  fx=p*dbinom(r-1,x-1,p)
  if(plotpdf==TRUE & plotcdf==FALSE)
  {
    discrete.plotpdf(x,fx)
  }
  if(plotpdf==FALSE & plotcdf==TRUE)
  {
    discrete.plotcdf(x,fx)
  }
  if(plotpdf==TRUE & plotcdf==TRUE)
  {
    par(mfrow=c(2,1))
    discrete.plotpdf(x,fx)
    discrete.plotcdf(x,fx)
  }
  par(mfrow=c(1,1))
  return(list(mean=mu,variance=sigma2,standard.deviation=sigma))
}

negbinom.prob=function(r,p,lb,ub=lb,inclusive="both")
{
  if(ub<lb){return("ub cannot be smaller than lb!")}
  if(ub<r){return(0)}
  if((as.integer(r)-r)!=0 | r<1){return("r must be an integer no less than 1")}
  if(p<=0 | p>1){return("p must be 0<p<1")}

  if(lb<r){a1=r}else{a1=floor(lb)}
  if(ub<Inf){
    b1=ceiling(ub)
    x=a1:b1
    fx=p*dbinom(r-1,x-1,p)
    if(inclusive=="both"){
      return(sum(fx[lb<=x &x<=ub]))
    }else if(inclusive=="left"){
      return(sum(fx[lb<=x &x<ub]))
    }else if(inclusive=="right"){
      return(sum(fx[lb<x &x<=ub]))
    }else if(inclusive=="none"){
      return(sum(fx[lb<x &x<ub]))
    }else{
      return("inclusive must be one of these: both, left, right, none")
    }
  }else if(ub==Inf){
    if(inclusive=="left" | inclusive=="both"){
      if(ceiling(lb)-1<r){return(1)
      }else{
          x=r:(ceiling(lb)-1)
          fx=p*dbinom(r-1,x-1,p)
          return(1-sum(fx))
        }
    }else if(inclusive=="right" | inclusive=="none"){
      if((floor(lb)<r)){return(1)
      }else{
          x=r:floor(lb)
          fx=p*dbinom(r-1,x-1,p)
          return(1-sum(fx))
        }
    }else{
      return("inclusive must be one of these: both, left, right, none")
    }
  }
}

###################################
#Geometric distribution
geometric.summary=function(p,plotpdf=TRUE,plotcdf=TRUE)
{
  if(p<=0 | p>1){return("p must be 0<p<1")}
  r=1
  mu=r/p
  sigma2=r*(1-p)/(p^2)
  sigma=sqrt(sigma2)

  x=r:max(r+10,ceiling(mu)+10)
  fx=p*dbinom(r-1,x-1,p)
  if(plotpdf==TRUE & plotcdf==FALSE)
  {
    discrete.plotpdf(x,fx)
  }
  if(plotpdf==FALSE & plotcdf==TRUE)
  {
    discrete.plotcdf(x,fx)
  }
  if(plotpdf==TRUE & plotcdf==TRUE)
  {
    par(mfrow=c(2,1))
    discrete.plotpdf(x,fx)
    discrete.plotcdf(x,fx)
  }
  par(mfrow=c(1,1))
  return(list(mean=mu,variance=sigma2,standard.deviation=sigma))
}

geometric.prob=function(p,lb,ub=lb,inclusive="both")
{
  if(ub<lb){return("ub cannot be smaller than lb!")}
  if(ub<1){return(0)}
  if(p<=0 | p>1){return("p must be 0<p<1")}
  r=1
  if(lb<r){a1=r}else{a1=floor(lb)}
  if(ub<Inf){
    b1=ceiling(ub)
    x=a1:b1
    fx=p*dbinom(r-1,x-1,p)
    if(inclusive=="both"){
      return(sum(fx[lb<=x &x<=ub]))
    }else if(inclusive=="left"){
      return(sum(fx[lb<=x &x<ub]))
    }else if(inclusive=="right"){
      return(sum(fx[lb<x &x<=ub]))
    }else if(inclusive=="none"){
      return(sum(fx[lb<x &x<ub]))
    }else{
      return("inclusive must be one of these: both, left, right, none")
    }
  }else if(ub==Inf){
    if(inclusive=="left" | inclusive=="both"){
      if(ceiling(lb)-1<r){return(1)
      }else{
        x=r:(ceiling(lb)-1)
        fx=p*dbinom(r-1,x-1,p)
        return(1-sum(fx))
      }
    }else if(inclusive=="right" | inclusive=="none"){
      if((floor(lb)<r)){return(1)
      }else{
        x=r:floor(lb)
        fx=p*dbinom(r-1,x-1,p)
        return(1-sum(fx))
      }
    }else{
      return("inclusive must be one of these: both, left, right, none")
    }
  }
}


###################################
#Hypergeometric distribution
hypergeo.summary=function(N,K,n,plotpdf=TRUE,plotcdf=TRUE)
{
  if((as.integer(N)-N)!=0 | N<=0){return("N must be a positive integer")}
  if((as.integer(K)-K)!=0 | K<=0){return("K must be a positive integer")}
  if((as.integer(n)-n)!=0 | n<=0){return("n must be a positive integer")}
  if(N<K){return("N must be >= K")}
  if(N<n){return("N must be >= n")}
  x=(max(0,n+K-N)):(min(K,n))
  fx=c()
  for(m in 1:length(x))
  {
    fx=c(fx,nCr(K,x[m])*nCr(N-K,n-x[m])/nCr(N,n))
  }

  mu=n*K/N
  sigma2=n*K/N*(1-K/N)*(N-n)/(N-1)
  sigma=sqrt(sigma2)

  if(plotpdf==TRUE & plotcdf==FALSE)
  {
    discrete.plotpdf(x,fx)
  }
  if(plotpdf==FALSE & plotcdf==TRUE)
  {
    discrete.plotcdf(x,fx)
  }
  if(plotpdf==TRUE & plotcdf==TRUE)
  {
    par(mfrow=c(2,1))
    discrete.plotpdf(x,fx)
    discrete.plotcdf(x,fx)
  }
  par(mfrow=c(1,1))
  return(list(mean=mu,variance=sigma2,standard.deviation=sigma))
}

hypergeo.prob=function(N,K,n,lb,ub=lb,inclusive="both")
{
  if(ub<lb){return("ub cannot be smaller than lb!")}
  if((as.integer(N)-N)!=0 | N<=0){return("N must be a positive integer")}
  if((as.integer(K)-K)!=0 | K<=0){return("K must be a positive integer")}
  if((as.integer(n)-n)!=0 | n<=0){return("n must be a positive integer")}
  if(N<K){return("N must be >= K")}
  if(N<n){return("N must be >= n")}
  x=(max(0,n+K-N)):(min(K,n))
  fx=c()
  for(m in 1:length(x))
  {
    fx=c(fx,nCr(K,x[m])*nCr(N-K,n-x[m])/nCr(N,n))
  }
  return(discrete.prob(x,fx,lb,ub,inclusive=inclusive))
}

###################################
#Poisson distribution
poisson.summary=function(lambda,L=1,plotpdf=TRUE,plotcdf=TRUE)
{
  if(lambda<=0){return("lambda must be >0")}
  if(L<=0){return("L must be >0")}
  mu=lambda*L
  sigma2=mu
  sigma=sqrt(sigma2)

  x=0:max(10,ceiling(mu)+10)
  fx=dpois(x,lambda*L)
  if(plotpdf==TRUE & plotcdf==FALSE)
  {
    discrete.plotpdf(x,fx)
  }
  if(plotpdf==FALSE & plotcdf==TRUE)
  {
    discrete.plotcdf(x,fx)
  }
  if(plotpdf==TRUE & plotcdf==TRUE)
  {
    par(mfrow=c(2,1))
    discrete.plotpdf(x,fx)
    discrete.plotcdf(x,fx)
  }
  par(mfrow=c(1,1))
  return(list(mean=mu,variance=sigma2,standard.deviation=sigma))
}

poisson.prob=function(lambda,L=1,lb,ub=lb,inclusive="both")
{
  if(ub<lb){return("ub cannot be smaller than lb!")}
  if(lambda<=0){return("lambda must be >0")}
  if(L<=0){return("L must be >0")}
  r=0
  if(lb<r){a1=r}else{a1=floor(lb)}
  if(ub<Inf){
    b1=ceiling(ub)
    x=a1:b1
    fx=dpois(x,lambda*L)
    if(inclusive=="both"){
      return(sum(fx[lb<=x &x<=ub]))
    }else if(inclusive=="left"){
      return(sum(fx[lb<=x &x<ub]))
    }else if(inclusive=="right"){
      return(sum(fx[lb<x &x<=ub]))
    }else if(inclusive=="none"){
      return(sum(fx[lb<x &x<ub]))
    }else{
      return("inclusive must be one of these: both, left, right, none")
    }
  }else if(ub==Inf){
    if(inclusive=="left" | inclusive=="both"){
      if(ceiling(lb)-1<r){return(1)
      }else{
        x=r:(ceiling(lb)-1)
        fx=dpois(x,lambda*L)
        return(1-sum(fx))
      }
    }else if(inclusive=="right" | inclusive=="none"){
      if((floor(lb)<r)){return(1)
      }else{
        x=r:floor(lb)
        fx=dpois(x,lambda*L)
        return(1-sum(fx))
      }
    }else{
      return("inclusive must be one of these: both, left, right, none")
    }
  }
}
