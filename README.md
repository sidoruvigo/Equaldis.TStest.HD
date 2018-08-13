# Equaldis.TStest.HD
A two-sample test for the equality of distributions for high-dimensional data

[SiDOR Group](http://sidor.uvigo.es/en/). [University of Vigo.](http://uvigo.gal/)

## Install the development version from GitHub
```r
devtools::install_github('sidoruvigo/Equaldis.TStest.HD')
```
## Description
This package implements four different tests proposed in Cousido-Rocha et al. (2018). These methods test the (global) null hypothesis of equality of the univariate marginals of the p-variate distributions in the two populations. In other words, the null hypothesis is an intersection of the p null hypothesis corresponding to p different two-sample problems. These methods are particularly well suited to the low sample size, high dimensional setting (n << p). The sample size can be as small as 2. The test accounts for the possibility that the p variables in each data set can be weakly dependent. Three of the methods arise from different approaches to estimate the variance of the same statistic. This statistic averages p individual statistics based on comparing the empirical characteristic functions computed from the two samples. The last method is an alternative global test whose statistic averages the p-values derived from applying permutation tests to the individual statistics mentioned above. When the global null hypothesis is rejected such permutation p-values can also be used to identify which variables have been contributed to this significance. The standarized version of each test statistic and its p-value are computed among other things.

## Authors
+ Cousido Rocha, Marta.
+ Soage González, José Carlos.
+ de Uña-Álvarez, Jacobo.
+ D. Hart, Jeffrey.

## References
Cousido-Rocha, M., de Uña-Álvarez J., and Hart, J. (2018). A two-sample test for the equality of distributions for high-dimensional data. Preprint.
