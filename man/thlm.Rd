% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/threshold-reg.R
\name{thlm}
\alias{thlm}
\title{Threshold regression with a censored covariate}
\usage{
thlm(formula, data, method = c("cc", "reverse", "deletion-threshold",
  "complete-threshold", "all"), B = 0, subset, x.upplim = NULL,
  t0 = NULL, control = thlm.control())
}
\arguments{
\item{formula}{A formula expression in the form \code{response ~ predictors}.
The response variable is assumed to be fully observed.
The \code{thlm} function can accommodate at most one censored covariate,
which is entered as an \code{Surv} object; see \code{survival::Surv} for more detail.
When all the covariates are uncensored, the \code{thlm} function returns a \code{lm} object.}

\item{data}{An optional data frame list or environment contains variables in the \code{formula} and
the \code{subset} argument. If left unspecified, the variables are taken
from \code{environment(formula)}, typically the environment from which \code{thlm} is called.}

\item{method}{A character string specifying the threshold regression methods to be used.
The following are permitted:
\describe{
  \item{\code{cc}}{for complete-cases regression}
  \item{\code{reverse}}{for reverse survival regression}
  \item{\code{deletion-threshold}}{for deletion threshold regression}
  \item{\code{complete-threshold}}{for complete threshold regression}
  \item{\code{all}}{for all four approaches}
}}

\item{B}{A numeric value specifies the bootstrap size for estimating
the standard deviation of regression coefficient for the censored
covariate when \code{method = "deletion-threshold"} or
\code{method = "complete-threshold"}.
When \code{B = 0}, only the beta estimate will be displayed.}

\item{subset}{An optional vector specifying a subset of observations to be
used in the fitting process.}

\item{x.upplim}{An optional numeric value specifies the upper support
of the censored covariate. When left unspecified, the maximum of the
censored covariate will be used.}

\item{t0}{An optional numeric value specifies the threshold when
\code{method = "dt"} or \code{"ct"}.
When left unspecified, an optimal threshold will be determined to optimize test power
using the proposed procedure in Qian et al (2018).}

\item{control}{A list of parameters. The parameters are
\describe{
  \item{\code{t0.interval}}{controls the end points of the interval to be searched for
the optimal threshold when \code{t0} is left unspecified}
  \item{\code{t0.plot}}{controls whether the objective function will be plotted.
When \code{t0.plot} is ture, both the raw \code{t0.plot} values and the smoothed estimates
(using local polynomial regression fitting) are plotted.}
}}
}
\description{
This function fits a linear regression model when there is
a censored covaraite. The method involves thresholding the continuous
covariate into a binary covariate. A collection of threshold
regression methods are implemented to obtain the estimator of the
regression coefficient as well as to test the significance of the
effect of the censored covariate. When there is no censoring,
the method reduces to the simple linear regression.

The model assumes the linear regression model:
\deqn{Y = a_0 + a_1X + a_2Z + e,} where X is the covariate of interest which
is subject to right censoring, Z is a covariate matrix that are fully
observed, Y is the response variable, and e is an independent randon
error term with mean 0 and finite variance.

The hypothesis test of association is based on the significance of the
regression coefficient, a1. However, when deletion threshold
regression or complete threshold regression is executed, an equivalent
but easy-to-evaluate test is performed. Namely, given a threshold
t*, we define a derived binary covariate, X*, such that X* = 1 when X
> t* and X* = 0 when X is uncensored and X < t*. The proposed linear
regression can be expressed as \deqn{E(Y|X^\ast, Z) = b_0 + b_1X^\ast +
b_2Z.} The proposed hypothesis test of association can be tested by the
significance of b1. Under the assumption that X is independent of Z
given X*, b2 is equivalent to a2.
}
\examples{
simDat <- function(n) {
  X <- rexp(n, 3)
  Z <- runif(n, 1, 6)
  Y <- 0.5 + 0.5 * X - 0.5 * Z + rnorm(n, 0, .75)
  cstime <- rexp(n, .75)
  delta <- (X <= cstime) * 1
  X <- pmin(X, cstime)
  data.frame(Y = Y, X = X, Z = Z, delta = delta)
}

set.seed(0)
dat <- simDat(200)

library(survival)
## Falsely assumes all covariates are free of censoring
thlm(Y ~ X + Z, data = dat)

## Complete cases regression
thlm(Y ~ Surv(X, delta) + Z, data = dat, method = "cc")

## reverse survival regression
thlm(Y ~ Surv(X, delta) + Z, data = dat, method = "rev")

## threshold regression without bootstrap
thlm(Y ~ Surv(X, delta) + Z, data = dat, method = "del")
thlm(Y ~ Surv(X, delta) + Z, data = dat, method = "com", control =
list(t0.interval = c(0.2, 0.6), t0.plot = FALSE))

## threshold regression with bootstrap
thlm(Y ~ Surv(X, delta) + Z, data = dat, method = "del", B = 100)
thlm(Y ~ Surv(X, delta) + Z, data = dat, method = "com", B = 100)

## display all
thlm(Y ~ Surv(X, delta) + Z, data = dat, method = "all", B = 100)
}
\references{
Qian, J., Chiou, S.H., Maye, J.E., Atem, F., Johnson, K.A. and Betensky, R.A. (2018)
Threshold regression to accommodate a censored covariate, \emph{Biometrics}, \bold{74}(4):
1261--1270.

Atem, F., Qian, J., Maye J.E., Johnson, K.A. and Betensky, R.A. (2017),
Linear regression with a randomly censored covariate: Application to an Alzheimer's study.
\emph{Journal of the Royal Statistical Society: Series C}, \bold{66}(2):313--328.
}
