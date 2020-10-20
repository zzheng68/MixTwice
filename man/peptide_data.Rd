\name{peptide_data}
\alias{peptide_data}
\docType{data}
\title{Peptide array data example}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("peptide_data")}
\format{
  A data frame with 152603 observations on the following 16 variables.
  \describe{
    the first 8 columns are RA patients and the rest columns are controls subjects
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
Zheng, Zihao, et al. \emph{Disordered Antigens and Epitope Overlap Between Anti Citrullinated Protein Antibodies and Rheumatoid Factor in Rheumatoid Arthritis.} \bold{Arthritis & Rheumatology} 72.2 (2020): 262-272.
}
\references{
\url{https://onlinelibrary.wiley.com/doi/abs/10.1002/art.41074}
}
\examples{
#### load the RA data

data(peptide_data)

#### pre-process the data

## each row is a peptide and each column is a subject; try rownames(peptide_data),colnames(peptide_data) 

## effect size

thetaHat = apply(peptide_data[1:8], 1, mean) - apply(peptide_data[9:16], 1, mean)

## standard error

sd1 = apply(peptide_data[,1:8], 1, sd)
sd0 = apply(peptide_data[,9:16], 1, sd)

s2 = sd1^2/8 + sd0^2/8

#### Apply on MixTwice

mm = mixtwice(thetaHat = thetaHat, s2 = s2, 
                  Btheta = 15, Bsigma = 10, df = 14, prop = 0.01)

plot(mm$grid.theta, cumsum(mm$mix.theta), type = "s", 
     xlab = "grid.theta", ylab = "ecdf of theta", lwd = 2)
plot(mm$grid.sigma, cumsum(mm$mix.sigma), type = "s",
     xlab = "grid.sigma", ylab = "ecdf of sigma", lwd = 2)

summary(mm$lfsr)

}
\keyword{datasets}