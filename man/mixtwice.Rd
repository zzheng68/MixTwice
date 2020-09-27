\name{mixtwice}
\alias{mixtwice}
\title{Large-scale hypothesis testing for high density peptide array by variance mixting}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mixtwice(theta, s2, Btheta, Bsigma, df, prop)
}
\arguments{
  \item{theta}{Estimated effect size}
  \item{s2}{Estimated standard error for effect size}
  \item{Btheta}{Half number of grid to estimated the distribution of effect size}
  \item{Bsigma}{Number of grid to estimated the distribution of standard error}
  \item{df}{Degree of freedom in chisquare distribution of estimated standard error}
  \item{prop}{Proportion of units randomly selected to fit the distribution}
}
\details{
\code{mixtwice} is a R function that takes estimated effect size and its standard error as input. It uses nonparametric MLE to fit the underlying distribution of effect size and standard error. By "mixing" on both effect size and variance, it could achieve promising statistical properties (i.e., FDR control, power, reproducibility, robustness) in large-scale testing.
}
\value{
\item{grid.theta}{Grid of distribution of effect size}
\item{grid.sigma}{Grid of distribution of standard error}
\item{prior.theta}{Estimated distribution of effect size, on grid.theta}
\item{prior.sigma}{Estimated distribution of standard error, on grid.sigma}
\item{lfdr}{Local false discovery rate for each testing unit}
\item{lfsr}{Local false sign rate for each testing unit}
}
\references{
cite our manuscript (on either github or ArXiv)
}
\author{
Zihao Zheng, Miriam A.Shelef, Michael A.Newton
}
\note{
further notes~~
}

\seealso{
See Also as \code{\link{help(array1)}}, \code{\link{help(array2)}}
}
\examples{
#### load the RA data

data(array1); data(array2)

#### pre-process the data

## effect size

delta1 = apply(array1[1:12], 1, mean) - apply(array1[13:24], 1, mean)
delta2 = apply(array2[1:8], 1, mean) - apply(array2[9:16], 1, mean)

## standard error

sd11 = apply(array1[,1:12], 1, sd)
sd10 = apply(array1[,13:24], 1, sd)

sd21 = apply(array2[,1:8], 1, sd)
sd20 = apply(array2[,9:16], 1, sd)

sigma1 = sqrt(sd11^2/12 + sd10^2/12)
sigma2 = sqrt(sd21^2/8 + sd20^2/8)

## z-score

get_zscore = function(x){
  
  x1 = x[1:length(x)/2]
  x2 = x[(length(x)/2+1):length(x)]
  
  t = t.test(x1, x2, var.equal = T)$statistic
  
  z = qnorm(pt(t, df = length(x)-2))
  
  return(z)
}

z1 = apply(array1, 1, get_zscore)
z1[z1==Inf] = max(z1[z1==Inf])

z2 = apply(array2, 1, get_zscore)

## p-value

p1 = 2*(1-pnorm(abs(z1)))
p2 = 2*(1-pnorm(abs(z2)))

#### Apply on MixTwice

library(alabama)

mixtwice1 = mixtwice(theta = delta1, s2 = sigma1^2, 
                 Btheta = 15, Bsigma = 10, df = 22, prop = 0.01)$lfsr

mixtwice2 = mixtwice(theta = delta2, s2 = sigma2^2, 
                 Btheta = 15, Bsigma = 10, df = 14, prop = 0.01)$lfsr

#### Apply on other large-scale testing

## ASH

library(ashr)

ash1 = get_lfsr(ash(delta1, sigma1, df = 22))

ash2 = get_lfsr(ash(delta2, sigma2, df = 14))

## Efron's local fdr

library(fdrtool)

lfdr1 = fdrtool(p1, statistic = "pvalue", plot = F)$lfdr

lfdr2 = fdrtool(p2, statistic = "pvalue", plot = F)$lfdr


## Storey's qvalue

library(qvalue)

qvalue1 = qvalue(p1)$lfdr
qvalue2 = qvalue(p2)$lfdr

## BH

BH1 = p.adjust(p1, method = "BH")
BH2 = p.adjust(p2, method = "BH")

## two-step ASH (Lu and Stephens 2019)
library(edgeR)
library(limma)

condition = rep(c("A","B"),each = 12 )
colnames(array1) = condition

dgecounts = calcNormFactors(DGEList(counts=array1, group=condition))

design = model.matrix(~condition)
v = voom(dgecounts,design,plot=FALSE)

lim = lmFit(v)
lim = eBayes(lim)
betahat = lim$coefficients[,2]
se = lim$stdev.unscaled[,2]*sqrt(lim$s2.post)
df = lim$df.total[1]

ebash1 = get_lfsr(ash(betahat, se, df=df))

condition = rep(c("A","B"),each = 8 ) # group indicator
colnames(array2) = condition

dgecounts = calcNormFactors(DGEList(counts=array2, group=condition))

design = model.matrix(~condition)
v = voom(dgecounts,design,plot=FALSE)

lim = lmFit(v)
lim = eBayes(lim)
betahat = lim$coefficients[,2]
se = lim$stdev.unscaled[,2]*sqrt(lim$s2.post)
df = lim$df.total[1]

ebash2 = get_lfsr(ash(betahat, se, df=df))

#### Evalute the reproducibility

d1 = data.frame(ash1, lfdr1,  mixtwice1, qvalue1, BH1, ebash1)
d2 = data.frame(ash2, lfdr2,  mixtwice2, qvalue2, BH2, ebash2)

## a summay reproducibility function
s = function(l1, l2, cut, plot = T){
  
  n1 = l1<=cut
  n2 = l2<=cut
  nboth = n1&n2
  neither = n1|n2

  return(list(n1 = sum(n1), n2 = sum(n2), nboth = sum(nboth), neither = sum(neither)))
}

C = 10^(-seq(3,5, by = 0.1))

Frac.ash <- Numb.ash <- Frac.2d <- Numb.2d <- Frac.efron <- Numb.efron <- 
  Frac.storey <- Numb.storey <- Frac.BH <- Numb.BH <- Frac.ebash <- Numb.ebash <- NULL

for (j in 1:length(C)) {
  
  a.ash = s(ash1, ash2, cut = C[j], plot = F)
  a.2d = s(mixtwice1, mixtwice2, cut = C[j], plot = F)
  a.efron = s(lfdr1, lfdr2, cut = C[j], plot = F)
  a.storey = s(qvalue1, qvalue2, cut = C[j], plot = F)
  a.BH = s(BH1, BH2, cut = C[j], plot = F)
  a.ebash = s(ebash1, ebash2, cut = C[j], plot = F)
  
  Frac.ash[j] = a.ash$nboth/a.ash$neither
  Numb.ash[j] = a.ash$nboth
  
  Frac.2d[j] = a.2d$nboth/a.2d$neither
  Numb.2d[j] = a.2d$nboth
  
  Frac.efron[j] = a.efron$nboth/a.efron$neither
  Numb.efron[j] = a.efron$nboth
  
  Frac.storey[j] = a.storey$nboth/a.storey$neither
  Numb.storey[j] = a.storey$nboth
  
  Frac.BH[j] = a.BH$nboth/a.BH$neither
  Numb.BH[j] = a.BH$nboth
  
  Frac.ebash[j] = a.ebash$nboth/a.ebash$neither
  Numb.ebash[j] = a.ebash$nboth
  
}

plot(C, Numb.ash, log = "x", ylim = c(0,max(Numb.2d)), 
     cex = 1, col = "blue", lwd = 2,
     xlab = "Significance Level", ylab = "Number of Significance (In Common)")

lines(C, Numb.ash, lwd = 2, col = "blue")

points(C, Numb.2d, 
       cex = 1, col = "green")
lines(C, Numb.2d, lwd = 2, col = "green")

points(C, Numb.efron, 
       cex = 1, col = "orange")
lines(C, Numb.efron, lwd = 2, col = "orange")

points(C, Numb.storey, 
       cex = 1, col = "purple")
lines(C, Numb.storey, lwd = 2, col = "purple")

points(C, Numb.BH, 
       cex = 1, col = "black")
lines(C, Numb.BH, lwd = 2, col = "black")

points(C, Numb.ebash, 
       cex = 1, col = "pink")
lines(C, Numb.ebash, lwd = 2, col = "pink")

legend("topleft",
       legend = c("BH", "Storey", "LocFDR", "ASH", "MixTwice", "VL+EB+ASH"),
       col = c('black', "purple", "orange", "blue", "green", "pink"),
       lwd = 1)

plot(C, Frac.ash, log = "x", ylim = c(0,0.7), 
     cex = 1, col = "blue", lwd = 2,
     xlab = "Significance Level", ylab = "Common Fraction in Both Array")

lines(C, Frac.ash, lwd = 2, col = "blue")

points(C, Frac.2d, 
       cex = 1, col = "green")
lines(C, Frac.2d, lwd = 2, col = "green")

points(C, Frac.efron, 
       cex = 1, col = "orange")
lines(C, Frac.efron, lwd = 2, col = "orange")

points(C, Frac.storey, 
       cex = 1, col = "purple")
lines(C, Frac.storey, lwd = 2, col = "purple")

points(C, Frac.BH, 
       cex = 1, col = "black")
lines(C, Frac.BH, lwd = 2, col = "black")

points(C, Frac.ebash, 
       cex = 1, col = "pink")
lines(C, Frac.ebash, lwd = 2, col = "pink")

legend("topleft",
       legend = c("BH", "Storey", "LocFDR", "ASH", "MixTwice", "VL+EB+ASH"),
       col = c('black', "purple", "orange", "blue", "green", "pink"),
       lwd = 1)
}
\keyword{High density peptide array}
\keyword{Large-scale hypothesis testing}
\keyword{Mixing distribution}
\keyword{Local false discovery rate}
