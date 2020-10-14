\name{mixtwice}
\alias{mixtwice}
\title{Large-scale hypothesis testing by variance mixing}
\description{
MixTwice deploys large-scale hypothesis testing in the case when testing units provide estimated effects and
estimated standard errors.  It produces empirical Bayesian local false discovery and sign rates for tests of zero effect.
}
\usage{
mixtwice(thetaHat, s2, Btheta, Bsigma, df, prop)
}
\arguments{
  \item{thetaHat}{Estimated effect size (vector over testing units))}
  \item{s2}{Estimated standard errors of thetaHat}
  \item{Btheta}{Grid size parameter for effect distribution}
  \item{Bsigma}{Grid size parameter for variance distribution}
  \item{df}{Degrees of freedom in chisquare associated with  estimated standard error}
  \item{prop}{Proportion of units randomly selected to fit the distribution}
}
\details{
\code{mixtwice} is a R function that takes estimated effects and standard errors. To compute local error-rate statistics, it finds nonparametric  MLEs of the underlying distributioons.
}
\value{
\item{grid.theta}{Support of the estimated mixing distribution on effects}
\item{grid.sigma}{Support of the estimated mixing distribution of variances}
\item{mix.theta}{Estimated distribution of effect size, on grid.theta}
\item{mix.sigma}{Estimated distribution of standard error, on grid.sigma}
\item{lfdr}{Local false discovery rate for each testing unit}
\item{lfsr}{Local false sign rate for each testing unit}
}
\references{
cite our manuscript (on either github or ArXiv)
}
\author{
Zihao Zheng,  Michael A.Newton
}
\note{
cite the biorxiv/arxiv paper
}

\seealso{
See Also as \code{?peptide_data}
}
\examples{
#### load the RA data

data(array)

#### pre-process the data

## effect size

delta = apply(array[1:8], 1, mean) - apply(array[9:16], 1, mean)

## standard error

sd1 = apply(array[,1:8], 1, sd)
sd0 = apply(array[,9:16], 1, sd)

sigma = sqrt(sd1^2/8 + sd0^2/8)

#### Apply on MixTwice

library(alabama)

result = mixtwice(thetaHat = delta, s2 = sigma^2, 
                 Btheta = 15, Bsigma = 10, df = 14, prop = 0.01)
}
\keyword{High density peptide array}
\keyword{Large-scale hypothesis testing}
\keyword{Mixing distribution}
\keyword{Local false discovery rate}
