\name{Discrete Distributions}
\alias{discrete.plotpdf}
\alias{discrete.plotcdf}
\alias{discrete.summary}
\alias{duniform.summary}
\alias{duniform.prob}
\alias{binomial.summary}
\alias{binomial.prob}
\alias{geometric.summary}
\alias{geometric.prob}
\alias{negbinom.summary}
\alias{negbinom.prob}
\alias{hypergeo.summary}
\alias{hypergeo.prob}
\alias{poisson.summary}
\alias{poisson.prob}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Discrete distributions
}
\description{
Plot pmf/cdf,  calculate mean/variance/standard deviation, and compute probabilities for various discrete distributions
}
\usage{
# For use-defined discrete distribution
discrete.plotpdf(x,fx)
discrete.plotcdf(x,fx)
discrete.summary(x,fx,plotpdf=c("TRUE","FALSE"), plotcdf=c("TRUE","FALSE"))
discrete.prob(x,fx,lb)
discrete.prob(x,fx,lb,ub,inclusive=c("none","left","right","both"))

# Discrete Uniform Distribution
duniform.summary(range,plotpdf=c("TRUE","FALSE"), plotcdf=c("TRUE","FALSE")))
duniform.prob(range,lb)
duniform.prob(range,lb,ub,inclusive=c("none","left","right","both"))

# Binomial distribution
binomial.summary(n,p,plotpdf=c("TRUE","FALSE"), plotcdf=c("TRUE","FALSE")))
binomial.prob(n,p,lb)
binomial.prob(n,p,lb,ub,inclusive=c("none","left","right","both"))

# Geometric distribution
geometric.summary(p,plotpdf=c("TRUE","FALSE"), plotcdf=c("TRUE","FALSE")))
geometric.prob(p,lb)
geometric.prob(p,lb,ub,inclusive=c("none","left","right","both"))

# Negative Binomial distribution
negbinom.summary(r,p,plotpdf=c("TRUE","FALSE"), plotcdf=c("TRUE","FALSE")))
negbinom.prob(r,p,lb)
negbinom.prob(r,p,lb,ub,inclusive=c("none","left","right","both"))

# Hypergeometric distribution
hypergeo.summary(N,K,n,plotpdf=c("TRUE","FALSE"), plotcdf=c("TRUE","FALSE")))
hypergeo.prob(N,K,n,lb)
hypergeo.prob(N,K,n,lb,ub,inclusive=c("none","left","right","both"))

# Poisson distribution
poisson.summary(lambda,L,plotpdf=c("TRUE","FALSE"), plotcdf=c("TRUE","FALSE")))
poisson.prob(lambda,L,lb)
poisson.prob(lambda,L,lb,ub,inclusive=c("none","left","right","both"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{possible values of a user defined discrete random variable}
  \item{fx}{probabilities of X=x, the order of entries in fx must matches the order of entries in x.}
  \item{plotpdf}{TRUE or FALSE, if TRUE, it plots the pmf}
  \item{plotcdf}{TRUE or FALSE, if TRUE, it polts the cdf}
  \item{lb,ub}{lower bound (lb) and upper bound (ub) in a probability statement; lb could be -Inf; ub could be Inf; ub cannot be less than lb}
  \item{inclusive}{"none": lb<X<ub;  "left": lb<=X<ub;  "right": lb<X<=ub;  "both": lb<=X<=ub}
  \item{Range}{contains all possible values of a discrete uniform random variable}
  \item{p}{parameter p of a Binomial distribution/Geometric distribution/Negtive Binomial distribution}
  \item{n}{parameter n of a Binomial distribution}
  \item{r}{parameter r of a negative Binomial distribution}
  \item{N,K,n}{parameters N, K, and n of a hypergeometric distribution}
  \item{lambda,L}{parameters lambda and L of a Poisson distribution. Default: L=1}
}
\details{
Plot the probability mass function (pmf) and the cumulative distribution function (cdf) and calculate probability, mean, variable and standard deviation, and compute probabilities of various discrete distributions (a self-defined discrete distribution, discrete uniform distribution, binomial distribution, geometric distribution, negative binomial distribution, hypergeometric distribution, and Poisson distribution).}
\value{
\item{discrete.plotpdf}{a figure of the pmf of the user-defined discrete distribution}
\item{discrete.plotcdf}{a figure of the cdf of the user-defined discrete distribution}
\item{discrete.summary}{a list of the mean, variance, and standard deviation of the user-defined discrete distribution, plot of the pmf/cdf or not}
\item{discrete.prob}{probability of X between lb and ub based on the user-defined discrete distribution}
\item{duniform.summary}{a list of the mean, variance, and standard deviation of a discrete uniform distribution, plot of the pmf/cdf or not}
\item{duniform.prob}{probability of X between lb and ub based on a uniform distribution}
\item{binomial.summary}{a list of the mean, variance, and standard deviation of a binomial distribution, plot of the pmf/cdf or not}
\item{binomial.prob}{probability of X between lb and ub based on a binomial distribution}
\item{geometric.summary}{a list of the mean, variance, and standard deviation of a geometric distribution, plot of the pmf/cdf or not}
\item{geometric.prob}{probability of X between lb and ub based on a geometric distribution}
\item{negbinom.summary}{a list of the mean, variance, and standard deviation of a negative binomial distribution, plot of the pmf/cdf or not}
\item{negbinom.prob}{probability of X between lb and ub based on a negative binomial distribution}
\item{hypergeo.summary}{a list of the mean, variance, and standard deviation of a hypergeometric distribution, plot of the pmf/cdf or not}
\item{hypergeo.prob}{probability of X between lb and ub based on a hypergeometric distribution}
\item{poisson.summary}{a list of the mean, variance, and standard deviation of a Poisson distribution, plot of the pmf/cdf or not}
\item{poisson.prob}{probability of X between lb and ub based on a Poisson distribution}
}

\references{
Chapter 3 of the textbook "Applied Statistics and Probability for Engineers" 7th edition
}
\author{
Dewei Wang
}
\note{
deweiwang@stat.sc.edu
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
x=c(0,1,2,3,4);fx=c(0.6561,0.2916,0.0486,0.0036,0.0001);

discrete.summary(x,fx)
discrete.summary(x,fx,plotpdf=FALSE,plotcdf=TRUE)
discrete.summary(x,fx,plotpdf=TRUE,plotcdf=FALSE)

discrete.prob(x,fx,2) #P(X=2)

discrete.prob(x,fx,2,4,inclusive="none") #P(2<X<4)
discrete.prob(x,fx,2,4,inclusive="left") #P(2<=X<4)
discrete.prob(x,fx,2,4,inclusive="right") #P(2<X<=4)
discrete.prob(x,fx,2,4,inclusive="both") #P(2<=X<=4)

discrete.prob(x,fx,2,Inf,inclusive="none") #P(2<X)
discrete.prob(x,fx,2,Inf,inclusive="left") #P(2<=X)
discrete.prob(x,fx,2,Inf,inclusive="right") #P(2<X)
discrete.prob(x,fx,2,Inf,inclusive="both") #P(2<=X)

discrete.prob(x,fx,-Inf,4,inclusive="none") #P(X<4)
discrete.prob(x,fx,-Inf,4,inclusive="left") #P(X<4)
discrete.prob(x,fx,-Inf,4,inclusive="right") #P(X<=4)
discrete.prob(x,fx,-Inf,4,inclusive="both") #P(X<=4)

#X~Discrete Uniform(1,2,4,5,8,10)
range=c(1,2,4,5,8,10)
duniform.summary(range)
duniform.prob(range,2) #P(X=2)
duniform.prob(range,2,4,inclusive="left") #P(2<=X<4)

#X~Binomial(n=5,p=0.3)
binomial.summary(5,0.3)
binomial.prob(5,0.3,2) #P(X=2)
binomial.prob(5,0.3,2,4,inclusive="left") #P(2<=X<4)

#X~Geometric(p=0.3)
geometric.summary(0.3)
geometric.prob(0.3,2) #P(X=2)
geometric.prob(0.3,2,4,inclusive="left") #P(2<=X<4)

#X~Negative Binomial(r=3,p=0.3)
negbinom.summary(3,0.3)
negbinom.prob(3,0.3,5) #P(X=5)
negbinom.prob(3,0.3,5,7,inclusive="left") #P(5<=X<7)

#X~Hypergeometric(N=300,K=100,n=4)
hypergeo.summary(300,100,4)
hypergeo.prob(300,100,4,2) #P(X=2)
hypergeo.prob(300,100,4,2,4,inclusive="left") #P(2<=X<4)

#X~Poisson(lambda=2.3,L=5)
poisson.summary(2.3,5)
poisson.prob(2.3,5,10) #P(X=10)
poisson.prob(2.3,5,1,Inf,inclusive="left") #P(1<=X)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
