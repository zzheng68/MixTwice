\name{array1}
\alias{array1}
\docType{data}
\title{Data 1 for CCP+RF+ vs Control peptide array.}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("array1")}
\format{
  A data frame with 152603 observations on the following 24 variables.
  \describe{
    First 12 columns are CCP+RF+ RA patients and the last columns are healthy controls
  }
}
\details{
The peptides in this data set is same as that in array2. However, the signals are evaluated independently are from different subjects.
}
\source{
Zheng, Zihao, et al. \emph{Disordered Antigens and Epitope Overlap Between Anti Citrullinated Protein Antibodies and Rheumatoid Factor in Rheumatoid Arthritis.} \bold{Arthritis & Rheumatology} 72.2 (2020): 262-272.
}
\references{
\url{https://onlinelibrary.wiley.com/doi/abs/10.1002/art.41074}
}
\seealso{
See Also as \code{\link{help(array2)}}, \code{\link{help(mixtwice)}}
}
\examples{
## load the RA data

data(array1)

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

## visualize z-score

hist(z1, probability = T, col = "blue", 100)
lines(density(rnorm(10^5)), lty = 2, lwd = 2)
}
\keyword{datasets}
