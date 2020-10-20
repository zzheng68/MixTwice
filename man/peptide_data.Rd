\name{peptide_data}
\alias{peptide_data}
\docType{data}
\title{Peptide array data example}
\description{
A high-density peptide microarray example to identify peptides for which antibody binding levels differ between control subjects and rheumatoid arthritis (RA) patients expressing a specific disease marker combination (i.e., CCP+RF+ RA). 
}
\usage{data("peptide_data")}
\format{
  A data frame with 152603 observations on the following 16 variables.
  \describe{
    The first 8 columns are RA patients and the remaining columns are from control subjects.
  }
}
\details{
Each row of the data (\code{rownames(peptide_data)}) is a probed length-12 peptide and each column of the data (\code{colnames(peptide_data)}) is a subject with distinct pseudo sample ID. The binding value is doubly-log transformed using natural base to stabilize variance.
}
\source{
Zheng, Zihao, et al. \emph{Disordered Antigens and Epitope Overlap Between Anti Citrullinated Protein Antibodies and Rheumatoid Factor in Rheumatoid Arthritis.} \bold{Arthritis & Rheumatology} 72.2 (2020): 262-272.
}
\references{
Zheng et al. \emph{MixTwice: Large scale hypothesis testing for peptide arrays by variance mixing}. Technical Report, October 2020.

\url{https://onlinelibrary.wiley.com/doi/abs/10.1002/art.41074}
}
\examples{
#### load the RA data

data(peptide_data)

#### pre-process the data

## effect size

thetaHat = apply(peptide_data[1:8], 1, mean) - apply(peptide_data[9:16], 1, mean)

## standard error

sd1 = apply(peptide_data[,1:8], 1, sd)
sd0 = apply(peptide_data[,9:16], 1, sd)

s2 = sd1^2/8 + sd0^2/8

#### Apply on MixTwice

mm = mixtwice(thetaHat = thetaHat, s2 = s2, 
                  Btheta = 15, Bsigma2 = 10, df = 14, prop = 0.01)

plot(mm$grid.theta, cumsum(mm$mix.theta), type = "s", 
     xlab = "grid.theta", ylab = "ecdf of theta", lwd = 2)
plot(mm$grid.sigma2, cumsum(mm$mix.sigma2), type = "s",
     xlab = "grid.sigma2", ylab = "ecdf of sigma2", lwd = 2)

summary(mm$lfsr)

}
\keyword{datasets}