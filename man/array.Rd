\name{array}
\alias{array}
\docType{data}
\title{Peptide array data example}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("array")}
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
## load the RA data

data(array)

## z-score

get_zscore = function(x){
  
  x1 = x[1:length(x)/2]
  x2 = x[(length(x)/2+1):length(x)]
  
  t = t.test(x1, x2, var.equal = T)$statistic
  
  z = qnorm(pt(t, df = length(x)-2))
  
  return(z)
}

z = apply(array, 1, get_zscore)

## visualize z-score

hist(z, probability = T, col = "blue", 100)
lines(density(rnorm(10^5)), lty = 2, lwd = 2)
}
\keyword{datasets}
