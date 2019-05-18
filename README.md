
**censCov**
-----------

[![Project Status: Concept – Minimal or no implementation has been done yet, or the repository is only intended to be a limited example, demo, or proof-of-concept.](https://www.repostatus.org/badges/latest/concept.svg)](https://www.repostatus.org/#concept) [![license](https://img.shields.io/badge/license-GPL%20(%3E=%203)-lightgrey.svg)](https://choosealicense.com/) [![Last-changedate](https://img.shields.io/badge/last%20change-2019--05--18-yellowgreen.svg)](/commits/master)

### Linear Regression with a Randomly Censored Covariate

***censCov*** consists of implementations of threshold regression approaches for linear regression models with a covariate subject to random censoring, including deletion threshold regression and completion threshold regression. Reverse survival regression, which flip the role of response variable and the covariate, is also considered.

### Installation

You can install and load **censCov** from CRAN using

``` r
install.packages("censCov")
library(censCov)
```

You can install `censCov` from github with:

``` r
## install.packages("devtools")
devtools::install_github("stc04003/censCov")
```

### References:

Qian, J., Chiou, S.H., Maye, J.E., Atem, F., Johnson, K.A. and Betensky, R.A. (2018). Threshold regression to accommodate a censored covariate, *Biometrics*, **74**(4): 1261--1270.

Atem, F., Qian, J., Maye J.E., Johnson, K.A. and Betensky, R.A. (2017), Linear regression with a randomly censored covariate: Application to an Alzheimer's study. *Journal of the Royal Statistical Society: Series C*, **66**(2):313--328.
