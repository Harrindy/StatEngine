\name{Linear Regression}
\alias{lm.est}
\alias{lm.coef.interval}
\alias{lm.coef.test}
\alias{lm.partialFtest}
\alias{lm.modelcheck}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Linear Regression
}
\description{
Inferences for linear regression: confidence intervals and hypthesis testing on regression coefficients, confidence interval on the mean response at a given x and prediction interval on a future observation at a given x, Partial F-test, model adequacy checking.
}
\usage{
# fit must be an lm object; e.g., fit=lm(y~x).
lm.est(fit)

# based on the fit, compute CIs and tests on betas
lm.coef.CI(fit,level=0.95)
# The default value of hypo.beta is zero
lm.coef.test(fit,alpha=0.05,H1="two",hypo.beta=?)

#Partial F-test,
#fit.H0 is the lm object using the null model
#fit.ALL is the lm object using the full model
lm.partialFtest(fit.H0=lmH0,fit.ALL=lmALL,alpha=0.05)

# Model Adequacy checking
lm.modelcheck(fit)
# fit must be an lm object. VIFs are also printed.
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{an lm object}
  \item{level}{the confidence level}
  \item{alpha}{the significance level}
  \item{hypo.beta}{the hypothesized beta value that about to be tested}
  \item{fit.H0, fit.ALL}{the null model and full model in the partial F-test}
}
\details{
Inferences for linear regression: confidence intervals and hypthesis testing on regression coefficients, confidence interval on the mean response at a given x and prediction interval on a future observation at a given x, Partial F-test, model adequacy checking.
}
\value{
\item{lm.est}{the least squares estimates}
\item{interval}{As long as the function has "interval", the outcome are confidence intervals.}
\item{test}{As long as the function has "test", it conduct the hypothesis testing.}
\item{lm.modelcheck}{residual analysis and VIFs}
}

\references{
Chapters 11-12 of the textbook "Applied Statistics and Probability for Engineers" 7th edition
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
Example1=read.csv("https://raw.githubusercontent.com/Harrindy/StatEngine/master/Data/HydrocarbonPurity.csv")
x=Example1$HydrocarbonLevels
y=Example1$Purity
fit=lm(y~x)
lm.est(fit)
summary(fit)
lm.coef.interval(fit,level=0.95)
predict.lm(fit,new=data.frame(x=1),interval="confidence")
predict.lm(fit,new=data.frame(x=1),interval="prediction")


Example2=read.csv("https://raw.githubusercontent.com/Harrindy/StatEngine/master/Data/WireBond.csv")
head(Example2,2)
y=Example2$PullStrength
x1=Example2$WireLength
x2=Example2$DieHeight
fit=lm(y~x1+x2)
summary(fit)
lm.coef.interval(fit,level=0.95)
lm.coef.test(fit,alpha=0.05,H1="two")
predict.lm(fit,new=data.frame(x1=8,x2=275),interval="confidence")
predict.lm(fit,new=data.frame(x1=8,x2=275),interval="prediction")
data.summary(fit$residuals)
vif(fit)
lm.modelcheck(fit)

x3=x1^2
x4=x2^2
lmH0=lm(y~x1+x2)
lmALL=lm(y~x1+x2+x3+x4)
lm.partialFtest(fit.H0=lmH0,fit.ALL=lmALL,alpha=0.05)
summary(lm(y~x1+x2+x3))
summary(lm(y~x1+x2+x4))
summary(lm(y~x1+x2+x3+x4))

Example3=read.csv("https://raw.githubusercontent.com/Harrindy/StatEngine/master/Data/AirplaneSidewallPanels.csv")
head(Example3,2)
y=Example3$cost
x=Example3$lotsize
plot(x,y)
fit=lm(y~x)
lines(x,fit$fitted.values)
x2=x^2
fit.quad=lm(y~x+x2)
lines(x,fit.quad$fitted.values,col="blue")
summary(fit)
summary(fit.quad)


Example4=read.csv("https://raw.githubusercontent.com/Harrindy/StatEngine/master/Data/SurfaceFinishData.csv")
head(Example4,2)
plot(Example4)
y=Example4$SurfaceFinish
x1=Example4$RPM
x2=Example4$TypeofCuttingTool
fit=lm(y~x1+x2)
summary(fit)
par(mar=c(4,4,.1,.1))
par(mfrow=c(1,2))
plot(x1[x2==0],y[x2==0])
lines(x1[x2==0],fit$fitted.values[x2==0])
plot(x1[x2==1],y[x2==1])
lines(x1[x2==1],fit$fitted.values[x2==1])

x3=x1*x2
fit2=lm(y~x1+x2+x3)
summary(fit2)
par(mfrow=c(1,1))
plot(x1,y)
lines(x1[x2==0],fit$fitted.values[x2==0])
lines(x1[x2==1],fit$fitted.values[x2==1])
lines(x1[x2==0],fit2$fitted.values[x2==0],col="blue")
lines(x1[x2==1],fit2$fitted.values[x2==1],col="blue")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

