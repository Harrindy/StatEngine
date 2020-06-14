\name{Descriptive Statistics}
\alias{DataSummary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Descriptive Statistics
}
\description{
Calculate and plot descriptive statistics
}
\usage{
data.summary(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the sample}
  \item{plot}{True of False; plot or not the stem and leaf diagram, scatterplot, histogram, boxplot, QQ-plot}
}
\details{
Calculate various numerical descpritive statistics and plot the stem and leaf diagram, scatterplot, histogram, boxplot, QQ-plot
}
\value{
\item{data.summary}{a table of numerical descpritive statistics (std stands for standard deviation), a stem and leaf diagram, and a figure consisting of scatterplot, histogram, boxplot, QQ-plot of the sample}
}

\references{
Chapter 6 of the textbook "Applied Statistics and Probability for Engineers" 7th edition
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
x=scan("https://raw.githubusercontent.com/Harrindy/StatEngine/master/Data/CompressiveStrength.csv")
data.summary(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
