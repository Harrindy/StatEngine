\name{One-Sample Confidence Intervals}
\alias{Zinterval}
\alias{sample.size.Zinterval}
\alias{Tinterval}
\alias{AZinterval}
\alias{Chi2interval}
\alias{Propinterval}
\alias{sample.size.Propinterval}
\alias{Predinterval}
\alias{t.quantile}
\alias{Chi2.quantile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
One-Sample Confidence Intervals
}
\description{
Compute confidence intervals on the population mean, population variance, and population proportaion. In addition, it computes prediction interval for a single future observation from a normal distribution.
}
\usage{
#CI for pupulation mean of a normal distribution when the population variance is known:
Zinterval(level,sigma,data) # if data available
Zinterval(level,sigma,n,barx) # if stats are provided
# Choice of sample size for estimating the population mean when error is specified
sample.size.Zinterval(level,sigma,E)

#CI for pupulation mean when the population variance is unknown and the distribution is normal
#or when the sample size is smaller than 25:
Tinterval(level,data) # if data available
Tinterval(level,n,barx,s) # if stats are provided

#Large-sample CI for pupulation mean:
AZinterval.data(level,data) # if data available
AZinterval.data(level,n,barx,s) # if stats are provided

#CI for pupulation variance (or standard deviation) of a normal distribution:
Chi2interval(level,data) # if data available
Chi2interval(level,n,s) # if stats are provided

#Large-sample CI for a pupulation proportion:
Propinterval(level,n,X)
# Choice of sample size for estimating a population proportion when error is specified
sample.size.Propinterval(level,ini.p,E) # using an intial guess
sample.size.Propinterval(level,ini.p=0.5,E) # using the conservative apporach

# Prediction interval of a single future observation form a normal distribution:
Predinterval(level,data) # if data available
Predinterval(level,n,barx,s) # if stats are provided
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{level}{the confidence level}
  \item{data}{a vector of the observed sample}
  \item{sigma}{the known population standard deviation}
  \item{s}{the observed sample standard deviation}
  \item{barx}{the observed sample mean}
  \item{n}{the sample size}
  \item{X}{number of observations belongs to a class of interest}
  \item{E}{specified error in sample size calculation}
  \item{ini.p}{A initial estimate of the populatin proportion. Default is 0.5 which corresponds to the conservative approach}
  \item{df}{the degrees of freedom of a t or chi.square distribution}
  \item{q}{a quantile value}
}
\details{
Compute CIs for the population mean, population variance, and population proportaion and PI for a single future observation from a normal distribution.}
\value{
\item{interval}{As long as the function has "interval", the outcome contains a two-sided CI (or PI) and the two one-sided confidence bounds.}
\item{sample.size}{As long as the function has "sample.size", the outcome is the minimum sample size required to control the error to be no larger than E.}
\item{t.quantile}{quantile of a t distribution}
\item{Chi2.quantile}{quantile of a chi.square distribution}
}

\references{
Chapter 8 of the textbook "Applied Statistics and Probability for Engineers" 7th edition
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
#Zinterval
# Useing data, must include the = sign
x=c(64.1, 64.7, 64.5, 64.6, 64.5, 64.3, 64.6, 64.8, 64.2, 64.3)
Zinterval(level=0.95,sigma=1,data=x)
Zinterval(level=0.99,sigma=1,data=x)
sample.size.Zinterval(E=0.5,sigma=1,level=0.95)
# Useing stats, must include the = sign
Zinterval(level=0.95,sigma=2,n=9,barx=98)

#Tinterval
# Useing stats, must include the = sign
Tinterval(level=0.95,n=10,barx=1000,s=20)
Tinterval(level=0.95,n=25,barx=1000,s=20)
Tinterval(level=0.99,n=10,barx=1000,s=20)
Tinterval(level=0.99,n=25,barx=1000,s=20)

#Large-sample Zinterval
# Useing data, must include the = sign
x=scan("https://raw.githubusercontent.com/Harrindy/StatEngine/master/Data/Mercury.csv")
AZinterval(level=0.95,data=x)

#Chi.square interval for variance/standard deviation
# Useing stats, must include the = sign
Chi2interval(level=0.95,n=20,s=0.01532)

#CIs for a porpulation proportion
#must include the = sign
Propinterval(level=0.95,n=85,X=10)
sample.size.Propinterval(level=0.95,ini.p=0.12,E=0.05)
sample.size.Propinterval(level=0.95,ini.p=0.5,E=0.05)

#Prediction interval for normal distribution
# Useing data, must include the = sign
x=c(19.8, 10.1, 14.9, 7.5, 15.4, 15.4, 15.4, 18.5, 7.9, 12.7, 11.9, 11.4, 11.4,
14.1, 17.6, 16.7, 15.8, 19.5, 8.8, 13.6, 11.9, 11.4)
Tinterval(level=0.95,data=x)
Predinterval(level=0.95,data=x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

