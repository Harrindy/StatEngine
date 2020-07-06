\name{One-Sample Confidence Intervals}
\alias{Ztest}
\alias{Ztest.power}
\alias{sample.size.Ztest}
\alias{Ttest}
\alias{Ttest.power}
\alias{sample.size.Ttest}
\alias{Chi2test}
\alias{Chi2test.power}
\alias{sample.size.Chi2test}
\alias{Proptest}
\alias{Proptest.power}
\alias{sample.size.Proptest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
One-Sample Tests
}
\description{
Conduct hypothesis testing on population mean, population variance, and population proportion. Compute power of a test, and calculate the required sample size for a desired power.
}
\usage{
#Test on pupulation mean of a normal distribution when the population variance is known:
Ztest(mu0,H1,alpha,sigma,sample) # if sample available
Ztest(mu0,h1,alpha,sigma,n,barx) # if statistics are provided
# Power calculation
Ztest.power(H1,alpha,sigma,n,delta)
# Sample size calculation
sample.size.Ztest(H1,sigma,alpha,beta,delta)

#Test on pupulation mean of a normal distribution when the population variance is unknown:
Ttest(mu0=?,H1=?,alpha=?,sample=?) #If data are available
Ttest(mu0=?,H1=?,alpha=?,n=?,barx=?,s=?) #If statistics are provided s=sd(data)
# Power calculation
Ttest.power(H1=?,est.sigma=sn,alpha=?,n=?,delta=?)
# Sample size calculation
sample.size.Ttest(H1=?,est.sigma=sn, beta=?,delta=?,alpha=?

#One-Sample Proportion Z-tests on a population proportion:
Proptest(p0=?,H1=?,alpha=?,n=?,X=?)
# Power calculation
Proptest.power(H1=?,alpha=?,p0=?,p1=?,n=?)
# Sample size calculation
sample.size.Proptest(H1=?,beta=?,alpha=?,p0=?,p1=?)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mu0}{the hypothesized value of the population mean in Z-tests and T-tests}
  \item{H1}{type of alternative hypothesis: "two", "left", or "right"}
  \item{alpha}{the significance level}
  \item{sample}{a vector of the observed sample}
  \item{sigma}{the known population standard deviation}
  \item{n}{the sample size}
  \item{barx}{the observed sample mean}
  \item{delta}{the delta value in power/sample size calculation of Z-test/T-test}
  \item{beta}{the desired power for sample size calculation}
  \item{s}{the observed sample standard deviation}
  \item{est.sigma}{an estimate of the population standard deviation in power/sample size calculation of T-test}
  \item{sigma0}{the hypothesized value of the population standard deviation in Chi-square-tests}
  \item{lambda}{the lambda value in power/sample size calculation of Chi-square-tests}
  \item{p0}{the hypothesized value of the population proportion}
  \item{X}{number of observations belongs to a class of interest}
  \item{p1}{the p1 value in power/sample size calculation of One-sample Prportion Z-tests}
}
\details{Conduct hypothesis testing on population mean, population variance, and population proportion. Compute power of a test, and calculate the required sample size for a desired power.}
\value{
\item{test}{As long as the function has "test", it produces the test results of using three approaches.}
\item{power}{As long as the function has "test", it computes the power of the test.}
\item{sample.size}{As long as the function has "sample.size", the outcome is the minimum sample size required to reach the given power.}
}

\references{
Chapter 9 of the textbook "Applied Statistics and Probability for Engineers" 7th edition
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
#Ztest
#must include the = sign
Ztest(mu0=3,H1="two",alpha=0.05,sigma=0.9,n=15,barx=2.78)
Ztest.power(H1="two",alpha=0.05,sigma=0.9,n=15,delta=3.25-3)
sample.size.Ztest(H1="two",sigma=0.9,alpha=0.05,beta=1-0.9,delta=3.75-3)

#Ttest
#must include the = sign
x=c(131.15, 130.69, 130.91, 129.54, 129.64, 128.77, 130.72, 128.33,
128.24, 129.65, 130.14, 129.29, 128.71, 129.00, 129.39, 130.42,
129.53, 130.12, 129.78, 130.92)
data.summary(x)
Ttest(mu0=130,H1="two",alpha=0.05,sample=x)
Ttest.power(H1="two",est.sigma=sd(x),alpha=0.05, n=length(x),delta=0.5)
sample.size.Ttest(H1="two",est.sigma =sd(x),beta=0.25,delta=0.1,alpha=0.05)

#Chi-square test
#must include the = sign
Chi2test(sigma0=sqrt(0.01),H1="right",alpha=0.05,n=20,s=sqrt(0.0153))
Chi2test.power(H1="right",alpha=0.05,n=20,lambda=1.25)
sample.size.Chi2test(H1="right",beta=0.2,lambda=1.25,alpha=0.05)

Chi2test(sigma0=0.01,H1="right",alpha=0.01,n=15,s=0.008)
Chi2test.power(H1="right",alpha=0.01,n=15,lambda=1.5)
sample.size.Chi2test(H1="right",beta=0.2,lambda=1.25,alpha=0.01)

#One-sample proportion Z-test
#must include the = sign
Proptest(p0=0.05, H1="left",alpha=0.05,n=300,X=7)
Proptest.power(H1="left",alpha=0.05,p0=0.05,p1=0.03,n=300)
sample.size.Proptest(H1="left",beta=0.1,alpha=0.05,p0=0.05,p1=0.03)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

