\name{Two-Sample Inferences}
\alias{twosample.Zinterval}
\alias{twosample.Ztest}
\alias{twosample.Tinterval}
\alias{twosample.Ttest}
\alias{Finterval}
\alias{Ftest}
\alias{twosample.Propinterval}
\alias{twosample.Proptest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Two-Sample Inferences (Confidence Intervals and Tests)
}
\description{
Compute confidence intervals and conduct hypothesis testing on difference between two population means, population variances, and population proportions.
}
\usage{
#CI for the difference in two pupulation means of a normal distribution
#when the population variances are known:
# if sample available
twosample.Zinterval(level,sigma1,sigma2,sample1,sample2)
# if stats are provided
twosample.Zinterval(level,sigma1,sigma2,barx1,barx2,n1,n2)

#Test on the difference in two pupulation means of a normal distribution
#when the population variances are known:
# if sample available
twosample.Ztest(Delta0,H1,alpha,sigma1,sigma2,sample1,sample2)
# if stats are provided
twosample.Ztest(Delta0,H1,alpha,sigma1,sigma2,barx1,barx2,n1,n2)

#CI for the difference in two pupulation means of a normal distribution
#when the population variances are unknown (pooled=yes or no)
# if sample available
twosample.Tinterval(level, pooled,sample1,sample2)
# if stats are provided
twosample.Tinterval(level, pooled,barx1,barx2,n1,n2,s1,s2)

#Test on the difference in two pupulation means of a normal distribution
#when the population variances are unknown:
# if sample available
twosample.Ttest(Delta0,H1,alpha,pooled=yes,sample1,sample2)
# if stats are provided
twosample.Ttest(Delta0,H1,alpha,pooled=yes,barx1,barx2,n1,n2,s1,s2)

#CI for the ratio between two pupulation variances of a normal distribution
# if sample available
Finterval(level, sample1,sample2)
# if stats are provided
Finterval(level, n1,n2,s1,s2)

#Test on the ratio between two pupulation variances of a normal distribution
# if sample available
Ftest(H1,alpha,sample1,sample2)
# if stats are provided
Ftest(H1,alpha,n1,n2,s1,s2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{level}{the confidence level}
  \item{sample1}{a vector of the observed sample from the first population}
  \item{sample2}{a vector of the observed sample from the second population}
  \item{sigma1}{the known population standard deviation of the first population}
  \item{sigma2}{the known population standard deviation of the second population}
  \item{s1,s2}{the sample standard deviations}
  \item{barx1,barx2}{the sample means}
  \item{n1,n2}{the sample sizes}
  \item{X1,X2}{number of observations belongs to a class of interest}
  \item{Delta0}{the hypothesized value of mu1-mu2}
  \item{H1}{type of alternative: "two","left", or "right"}
  \item{alpha}{the significance level}
  \item{pooled}{"yes" or "no"}
}
\details{
Compute confidence intervals and conduct hypothesis testing on difference between two population means, population variances, and population proportions.
}
\value{
\item{interval}{As long as the function has "interval", the outcome contains a two-sided CI and the two one-sided confidence bounds.}
\item{test}{As long as the function has "test", it produces the test results of using three approaches.}
}

\references{
Chapter 10 of the textbook "Applied Statistics and Probability for Engineers" 7th edition
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
#two-sample Zinterval
#must include the = sign
twosample.Zinterval(level=0.9,sigma1=1,sigma2=1.5,barx1=87.6,barx2=74.5,n1=10,n2=12)

#two-sample Ztest
twosample.Ztest(Delta0=0,H1="right",alpha=0.05, sigma1=8,sigma2=8,
                barx1=121,barx2=112,n1=10,n2=10)

#two-sample Tinterval
#must include the = sign
twosample.Tinterval(level=0.95,pooled="no",s1=5,s2=4,barx1=90,barx2=87,n1=10,n2=15)

#two-sample Ttest
catalyst1=c(91.50,94.18,92.18,95.39,91.79,89.07,94.72,89.21)
catalyst2=c(89.19,90.95,90.46,93.21,97.19,97.04,91.07,92.75)
data.summary(catalyst1)
data.summary(catalyst2)
twosample.Tinterval(level=0.95,pooled="yes",
                    sample1=catalyst1,sample2=catalyst2)
twosample.Ttest(Delta0=0,H1="two",alpha=0.05, pooled="yes",
                sample1=catalyst1,sample2=catalyst2)


C50=c(0.047, 0.060, 0.061, 0.064, 0.080, 0.090, 0.118, 0.165, 0.183)
C60=c(0.062, 0.105, 0.118, 0.137, 0.153, 0.197, 0.210, 0.250, 0.335)
data.summary(C50)
data.summary(C60)
twosample.Tinterval(level=0.95,pooled="no",
                    sample1=C50,sample2=C60)
twosample.Ttest(Delta0=0,H1="left",alpha=0.05, pooled="no",
                sample1=C50,sample2=C60)


#Paired T-test
Karlsrube=c(1.186,1.151,1.322,1.229,1.200,1.402,1.365,1.537,1.559)
   Lehigh=c(1.061,0.992,1.063,1.062,1.065,1.178,1.037,1.086,1.052)
   data.summary(Karlsrube-Lehigh)
Tinterval(level=0.95,sample=Karlsrube-Lehigh)
Ttest(mu0=0, H1="two",alpha=0.05,sample=Karlsrube-Lehigh)

#F-interval
Finterval(level=0.9,n1=11,n2=16,s1=5.1,s2=4.7)
Ftest(H1="two",alpha=0.1,n1=11,n2=16,s1=5.1,s2=4.7)

#two-sample proportion Z-interval
twosample.Propinterval(level=0.95, n1=85,n2=85,X1=10,X2=8)
twosample.Propinterval(level=0.95, n1=100,n2=100,X1=27,X2=19)
twosample.Proptest(H1="right",alpha=0.05,n1=100,n2=100,X1=27,X2=19)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

