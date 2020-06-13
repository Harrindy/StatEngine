\name{Continuous Distributions}
\alias{uniform.summary}
\alias{uniform.prob}
\alias{uniform.quantile}
\alias{normal.summary}
\alias{normal.prob}
\alias{normal.quantile}
\alias{exponential.summary}
\alias{exponential.prob}
\alias{exponential.quantile}
\alias{gamma.summary}
\alias{gamma.prob}
\alias{gamma.quantile}
\alias{weibull.summary}
\alias{weibull.prob}
\alias{weibull.quantile}
\alias{lognormal.summary}
\alias{lognormal.prob}
\alias{lognormal.quantile}
\alias{beta.summary}
\alias{beta.prob}
\alias{beta.quantile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Continuous distributions
}
\description{
Plot pdf/cdf,  calculate mean/variance/standard deviation, and compute probabilities for various continuous distributions
}
\usage{
# continuous Uniform Distribution
uniform.summary(a,b,plotpdf=c("TRUE","FALSE"), plotcdf=c("TRUE","FALSE")))
uniform.prob(a,b,lb,ub)
uniform.quantile(a,b,q)

# Normal distribution
normal.summary(mu,sigma,plotpdf=c("TRUE","FALSE"), plotcdf=c("TRUE","FALSE")))
normal.prob(mu,sigma,lb,ub)
normal.quantile(mu,sigma,q)

# Exponential distribution
exponential.summary(lambda,plotpdf=c("TRUE","FALSE"), plotcdf=c("TRUE","FALSE")))
exponential.prob(lambda,lb,ub)
exponential.quantile(lambda,q)

# Gamma distribution
gamma.summary(r,lambda,plotpdf=c("TRUE","FALSE"), plotcdf=c("TRUE","FALSE")))
gamma.prob(r,lambda,lb,ub)
gamma.quantile(r,lambda,p)

# Weibull distribution
weibull.summary(beta,delta,plotpdf=c("TRUE","FALSE"), plotcdf=c("TRUE","FALSE")))
weibull.prob(beta,delta,lb,ub)
weibull.quantile(beta,delta,p)

# Lognormal distribution
lognormal.summary(theta,omega,plotpdf=c("TRUE","FALSE"), plotcdf=c("TRUE","FALSE")))
lognormal.prob(theta,omega,lb,ub)
lognormal.quantile(theta,omega,p)

#Beta distribution
beta.summary(alpha,beta,plotpdf=c("TRUE","FALSE"), plotcdf=c("TRUE","FALSE")))
beta.prob(alpha,beta,lb,ub)
beta.quantile(alpha,beta,p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{plotpdf}{TRUE or FALSE, if TRUE, it plots the pmf}
  \item{plotcdf}{TRUE or FALSE, if TRUE, it polts the cdf}
  \item{lb,ub}{lower bound (lb) and upper bound (ub) in a probability statement; lb could be -Inf; ub could be Inf; ub cannot be less than lb}
  \item{a,b}{lower bound and upper bound of the support of a continuous uniform distribution}
  \item{mu,sigma}{mean and standard deviation of a normal distribution}
  \item{lambda}{parameter of an exponential distribution}
  \item{r,lambda}{shape and scale parameters of a gamma distribution}
  \item{beta,delta}{shape and scale parameters of a Weibull distribution}
  \item{theta,omega}{shape and scale parameters of a lognormal distribution}
  \item{alpba,beta}{parameters of a beta distribution}
}
\details{
Plot the probability density function (pdf) and the cumulative distribution function (cdf) and calculate probability, mean, variable and standard deviation, and compute probabilities of various discrete distributions (continuous uniform distribution, normal distribution, exponential distribution, gamma distribution, Weibull distribution, lognormal distribution, beta distribution).}
\value{
\item{uniform.summary}{a list of the mean, variance, and standard deviation of a continuous uniform distribution, plot of the pdf/cdf or not}
\item{uniform.prob}{probability of X between lb and ub based on a continuous uniform distribution}
\item{uniform.quantile}{quantile of a continuous uniform distribution}
\item{normal.summary}{a list of the mean, variance, and standard deviation of a normal distribution, plot of the pdf/cdf or not}
\item{normal.prob}{probability of X between lb and ub based on a normal distribution}
\item{normal.quantile}{quantile of a normal distribution}
\item{exponential.summary}{a list of the mean, variance, and standard deviation of an exponential distribution, plot of the pdf/cdf or not}
\item{exponential.prob}{probability of X between lb and ub based on an exponential distribution}
\item{exponential.quantile}{quantile of an exponential distribution}
\item{gamma.summary}{a list of the mean, variance, and standard deviation of a gamma distribution, plot of the pdf/cdf or not}
\item{gamma.prob}{probability of X between lb and ub based on a gamma distribution}
\item{gamma.quantile}{quantile of a gamma distribution}
\item{weibull.summary}{a list of the mean, variance, and standard deviation of a Weibull distribution, plot of the pdf/cdf or not}
\item{weibull.prob}{probability of X between lb and ub based on a Weibull distribution}
\item{weibull.quantile}{quantile of a Weibull distribution}
\item{lognormal.summary}{a list of the mean, variance, and standard deviation of a lognormal distribution, plot of the pdf/cdf or not}
\item{lognormal.prob}{probability of X between lb and ub based on a lognormal distribution}
\item{lognormal.quantile}{quantile of a lognormal distribution}
\item{beta.summary}{a list of the mean, variance, and standard deviation of a beta distribution, plot of the pdf/cdf or not}
\item{beta.prob}{probability of X between lb and ub based on a beta distribution}
\item{beta.quantile}{quantile of a beta distribution}
}

\references{
Chapter 4 of the textbook "Applied Statistics and Probability for Engineers" 7th edition
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
# Continuous uniform distribution
a=4.9;b=5.1;uniform.summary(a,b)
uniform.prob(a,b,4.95,5) # P(4.95<X<5)
uniform.quantile(a,b,0.9) # x such that P(X>x)=0.1

#Normal distribution
normal.summary(10,2) #mu=10,sigma=2,variance=4
normal.prob(10,2,9,11)
normal.quantile(10,2,0.98)

#Exponential distribution
exponential.summary(25) #lambda=25
exponential.prob(25,0.1,Inf) #P(X>0.1)
exponential.quantile(25,0.1) # x such that P(X>x)=0.9

#Gamma distribution
gamma.summary(10,.5) #r=10,lambda=0.5
gamma.prob(10,0.5,25,Inf) #P(X>25)
gamma.quantile(10,0.5,0.95) # x such that P(X<x)=0.95

#Weibull distribution
weibull.summary(2,5000) #beta=2,delta=5000
weibull.prob(2,5000,6000,Inf) #P(X>6000)
weibull.quantile(2,5000,1-0.05) # x such that P(X>x)=0.05

#Lognormal distribution
lognormal.summary(10,1.5) #theta=10,omega=1.5
lognormal.prob(10,1.5,10000,Inf) #P(X>10000)
lognormal.quantile(10,1.5,1-0.99) # x such that P(X>x)=0.99

#Beta distribution
beta.summary(2.5,1) #alpha=2.5, beta=1
beta.prob(2.5,1,0.7,Inf) # P(X>0.7)
beta.quantile(2.5,1,0.99) # x such that P(X<x)=0.99
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
