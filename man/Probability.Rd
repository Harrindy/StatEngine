\name{Counting Techniques}
\alias{nPr}
\alias{nCr}
\alias{SimPerm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sample Space Characterization
}
\description{
Use permutations and combinations to count outcomes
}
\usage{
nPr(n,r): permute r items from a set of n distinct items

nCr(n,r): select r items from a set of n items where order does not matter

SimPerm(n_vec): permute n items when they are not totally distinct
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{total n itmes}
  \item{r}{select r}
  \item{n_vec}{a vector of (n_1, n_2,...,n_r)}
}
\details{
Permutation, Combination, and Permutation of similar items
}
\value{
nPr(n,r) provides the number of different ways to permute r items from a set of n distinct items.

nCr(n,r) provides the number of different ways to select r items from a set of n items where order does not matter.

SimPerm(n_vec) provides the number of different ways to permute n items when they are not totally distinct.
}
\references{
Section 2.2 of the textbook "Applied Statistics and Probability for Engineers" 7th edition
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
nPr(8,4)
nCr(3,2)
nCr(47,4)
SimPerm(c(2,3,2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
