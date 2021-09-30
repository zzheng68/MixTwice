# MixTwice

This is a Github reporsitory for `MixTwice`

`MixTwice` is an empirical Bayes approach  for large-scale hypothesis testing.  It is currently developed for the two-group comparison problem, and takes as input a vector of estimated effects (unit-specific differences between groups, such as log fold change) and a second vector of estimated standard errors for these effects.   Using a shape-constrained, semi-parametric mixture model, it computes unit-specific local false discovery rates and local false sign rates, which may be used to prioritize units for follow-up analysis.   Read a more complete description at [Zheng et al, 2021, Bioinformatics](https://academic.oup.com/bioinformatics/article/37/17/2637/6162883).  

`MixTwice` is available in [CRAN](https://cran.rstudio.com/web/packages/MixTwice/index.html), and also here following the instructions below.


##### Install `MixTwice` R package: 

1. To locally download the `MixTwice` package, you can [use this link](https://github.com/wiscstatman/MixTwice/releases/tag/v1.1.1) to download the .zip file and install on R.

2. Install from Github

```R
install.packages("devtools")
devtools::install_github("wiscstatman/MixTwice")
```



