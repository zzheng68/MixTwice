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
  \item{thetaHat}{Estimated effect sizes (vector over testing units))}
  \item{s2}{Estimated squared standard errors of thetaHat}
  \item{Btheta}{Grid size parameter for effect distribution}
  \item{Bsigma2}{Grid size parameter for variance distribution}
  \item{df}{Degrees of freedom in chisquare associated with  estimated standard error}
  \item{prop}{Proportion of units randomly selected to fit the distribution, with default, \code{prop = 1} (use all units to fit the distribution).}
}
\details{
\code{mixtwice} takes estimated effects and standard errors over a set of testing units. To compute local error-rate statistics, it finds nonparametric  MLEs of the underlying distributions.
}
\value{
\item{grid.theta}{Support of the estimated mixing distribution on effects}
\item{grid.sigma2}{Support of the estimated mixing distribution of variances}
\item{mix.theta}{Estimated distribution of effect size, on grid.theta}
\item{mix.sigma2}{Estimated distribution of variance, on grid.sigma2}
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
### basic setting 

pi = 0.5 ## true value of null proportion

p = 1000; n = 10 ## number of testing units and sample size of each unit

p1=(1-pi)*p; p2=pi*p ## number of non-null and null

### Example 1: for a single group testing problem for zero effect

mu=c(rnorm(round(p1), mean=0, sd=1), rep(0, round(p2)))

sd=rep(1, p)

x=NULL

for (i in 1:(p1+p2)) {
  
  xx=rnorm(n, mu[i], sd[i])
  
  x=rbind(x,xx)
}


thetaHat = rowMeans(x) ## effect size of each testing unit

s2 = apply(x, 1, sd)^2/n ## estimated variance of effect size

mm1=mixtwice(thetaHat=thetaHat, s2=s2, Btheta = 15, Bsigma2 = 10, df=n-1)

## summarize and visualize the result

# estimated mixing distribution and true mixing distribution
plot(mm1$grid.theta, cumsum(mm1$mix.theta), type = "s", 
     xlab = "grid.theta", ylab = "ecdf of theta", lwd = 2)
lines(ecdf(mu),cex = 0.1, lwd = 0.5, lty = 2, col = "red")

plot(mm1$grid.sigma2, cumsum(mm1$mix.sigma2), type = "s",
     xlab = "grid.sigma2", ylab = "ecdf of sigma2", lwd = 2)

# effect size and estimated local false discovery rate

plot(thetaHat, mm1$lfdr, pch = ".", cex = 2, xlab = "estimated effect size", ylab = "local false discovery rate")

# true positive and false positive rate, under the level of 0.05

mean(mm1$lfsr[1:500]<=0.05) # true positive rate
mean(mm1$lfsr[501:1000]<=0.05) # false positive rate

# null proportion estimation

max(mm1$mix.theta)

### Example 2: for a two-group comparison testing problem for zero difference

mu0 = rep(0, p)
mu1=c(rnorm(round(p1), mean=0, sd=1), rep(0, round(p2)))

sd=rep(1, p)

x0<-x1<-NULL

for (i in 1:p) {
  
  xx0=rnorm(n, mu0[i], sd[i])
  
  x0=rbind(x0,xx0)
  
  xx1=rnorm(n, mu1[i], sd[i])
  
  x1=rbind(x1,xx1)
}


thetaHat = rowMeans(x1) - rowMeans(x0) ## effect size of each testing unit

s2 = apply(x1, 1, sd)^2/n + apply(x0, 1, sd)^2/n ## estimated variance of effect size

mm2=mixtwice(thetaHat=thetaHat, s2=s2, Btheta = 15, Bsigma2 = 10, df=2*n-2)

## summarize and visualize the result

# estimated mixing distribution and true mixing distribution
plot(mm2$grid.theta, cumsum(mm2$mix.theta), type = "s", 
     xlab = "grid.theta", ylab = "ecdf of theta", lwd = 2)
lines(ecdf(mu),cex = 0.1, lwd = 0.5, lty = 2, col = "red")

plot(mm2$grid.sigma2, cumsum(mm2$mix.sigma2), type = "s",
     xlab = "grid.sigma2", ylab = "ecdf of sigma2", lwd = 2)

# effect size and estimated local false discovery rate

plot(thetaHat, mm2$lfdr, pch = ".", cex = 2, xlab = "estimated effect size", ylab = "local false discovery rate")

# true positive and false positive rate, under the level of 0.05

mean(mm2$lfsr[1:500]<=0.05) # true positive rate
mean(mm2$lfsr[501:1000]<=0.05) # false positive rate

# null proportion estimation

pi_est = max(mm2$mix.theta)
}
\keyword{High density peptide array}
\keyword{Large-scale hypothesis testing}
\keyword{Mixing distribution}
\keyword{Local false discovery rate}
