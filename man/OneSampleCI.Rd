\name{One-Sample Inference}
\alias{Zinterval.data}
\alias{Zinterval.stat}
\alias{sample.size.Zinterval}
\alias{Tinterval.data}
\alias{Tinterval.stat}
\alias{t.quantile}
\alias{AZinterval.data}
\alias{AZinterval.stat}
\alias{Chi2interval.data}
\alias{Chi2interval.stat}
\alias{Chi2.quantile}
\alias{Propinterval}
\alias{sample.size.Propinterval}
\alias{Predinterval.data}
\alias{Predinterval.stat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
One-Sample Confidence Intervals and Hypothesis Tests
}
\description{
Compute confidence intervals and conduct hypothesis tests and on the population mean, population variance, and population proportaion. In addition, it computes prediction interval for a single future observation from a normal distribution.
}
\usage{
#CI for pupulation mean of a normal distribution when the population variance is known:
Zinterval.data(data,sigma,level)
Zinterval.data(n,barx,sigma,level)
# Choice of sample size for estimating the population mean when error is specified
sample.size.Zinterval(E,sigma,level)

#CI for pupulation mean when the population variance is unknown and the distribution is normal
#or when the sample size is smaller than 25:
Tinterval.data(data,level)
Tinterval.data(n,barx,s,level)

#Large-sample CI for pupulation mean:
AZinterval.data(data,level)
AZinterval.data(n,barx,s,level)

#CI for pupulation variance (or standard deviation) of a normal distribution:
Chi2interval.data(data,level)
Chi2interval.data(n,s,level)

#Large-sample CI for a pupulation proportion:
Propinterval(n,x,level)
# Choice of sample size for estimating a population proportion when error is specified
sample.size.Propinterval(E,ini.p,level)
sample.size.Propinterval(E,ini.p=0.5,level) # using the conservative apporach

# Prediction interval of a single future observation form a normal distribution:
Predinterval.data(data,level)
Predinterval.stat(n,barx,s,level)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{level}{the confidence level}
  \item{data}{a vector of the observed sample}
  \item{sigma}{the known population standard deviation}
  \item{s}{the observed sample standard deviation}
  \item{barx}{the observed sample mean}
  \item{n}{the sample size}
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
data=c(64.1, 64.7, 64.5, 64.6, 64.5, 64.3, 64.6, 64.8, 64.2, 64.3)
Zinterval.data(data,sigma=1,level=0.95)
Zinterval.data(data,sigma=1,level=0.99)
sample.size.Zinterval(E=0.5,sigma=1,level=0.95)
Zinterval.stat(n=9,barx=98,sigma=2,level=0.95)

Tinterval.stat(n=10,barx=1000,s=20,level=0.95)
Tinterval.stat(n=25,barx=1000,s=20,level=0.95)
Tinterval.stat(n=10,barx=1000,s=20,level=0.99)
Tinterval.stat(n=25,barx=1000,s=20,level=0.99)

data=scan("https://raw.githubusercontent.com/Harrindy/StatEngine/master/Data/Mercury.csv")
AZinterval.data(data,level=0.95)

Chi2interval.stat(n=20,s=0.01532,level=0.95)

Propinterval(n=85,x=10,level=0.95)
sample.size.Propinterval(E=0.05,ini.p=0.12,level=0.95)
sample.size.Propinterval(E=0.05,ini.p=0.5,level=0.95)

data=c(19.8, 10.1, 14.9, 7.5, 15.4, 15.4, 15.4, 18.5, 7.9, 12.7, 11.9, 11.4, 11.4,
14.1, 17.6, 16.7, 15.8, 19.5, 8.8, 13.6, 11.9, 11.4)
Tinterval.data(data,level=0.95)
Predinterval.data(data,level=0.95)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

