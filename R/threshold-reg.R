#' @name thlm
#' @rdname thlm
#' @title Threshold regression with a censored covariate
#'
#' @description This function fits a linear regression model when there is
#' a censored covaraite. The method involves thresholding the continuous
#' covariate into a binary covariate. A collection of threshold
#' regression methods are implemented to obtain the estimator of the
#' regression coefficient as well as to test the significance of the
#' effect of the censored covariate. When there is no censoring,
#' the method reduces to the simple linear regression.
#'
#' The model assumes the linear regression model:
#' \deqn{Y = a_0 + a_1X + a_2Z + e,} where X is the covariate of interest which
#' is subject to right censoring, Z is a covariate matrix that are fully
#' observed, Y is the response variable, and e is an independent randon
#' error term with mean 0 and finite variance.
#'
#' The hypothesis test of association is based on the significance of the
#' regression coefficient, a1. However, when deletion threshold
#' regression or complete threshold regression is executed, an equivalent
#' but easy-to-evaluate test is performed. Namely, given a threshold
#' t*, we define a derived binary covariate, X*, such that X* = 1 when X
#' > t* and X* = 0 when X is uncensored and X < t*. The proposed linear
#' regression can be expressed as \deqn{E(Y|X^\ast, Z) = b_0 + b_1X^\ast +
#' b_2Z.} The proposed hypothesis test of association can be tested by the
#' significance of b1. Under the assumption that X is independent of Z
#' given X*, b2 is equivalent to a2.
#' 
#' @param formula A formula expression in the form \code{response ~ predictors}.
#' The response variable is assumed to be fully observed.
#' See the documentation of \code{lm} or \code{formula} for more details.
#' @param data An optional data frame list or environment contains variables in the \code{formula} and
#' the \code{subset} argument. If left unspecified, the variables are taken
#' from \code{environment(formula)}, typically the environment from which \code{thlm} is called.
## #' @param cens An optional vector of censoring indicator (0 = censored, 1 = uncensored)
## #' for the censored covariate, which is assumed to be the first covariate specified in the
## #' \code{formula}. When \code{cens} is left unspecified or a vector of 1's,
## #' the model assumes all covariates are fully observed and the model reduced to simple
#' linear regression, i.e. \code{lm}.
#' @param subset An optional vector specifying a subset of observations to be
#' used in the fitting process.
#' @param method A character string specifying the threshold regression methods to be used.
#' The following are permitted:
#' \describe{
#'   \item{\code{cc}}{for complete-cases regression}
#'   \item{\code{rev}}{for reverse survival regression}
#'   \item{\code{dt}}{for deletion threshold regression}
#'   \item{\code{ct}}{for complete threshold regression}
#'   \item{\code{all}}{for all four approaches}
#' }
#' @param B A numeric value specifies the bootstrap size for estimating
#' the standard deviation of regression coefficient for the censored
#' covariate when \code{method = "dt"} or \code{"ct"}.
#' When \code{B = 0}, only the beta estimate will be displayed.
#' @param x.upplim An optional numeric value specifies the upper support
#' of censored covariate. When left unspecified, the maximum of the
#' censored covariate will be used.
#' @param t0 An optional numeric value specifies the threshold when
#' \code{method = "dt"} or \code{"ct"}.
#' When left unspecified, an optimal threshold will be determined to optimize test power
#' using the proposed procedure in Qian et al (2018).
#' @param control A list of parameters. The parameters are
#' \describe{
#'   \item{\code{t0.interval}}{controls the end points of the interval to be searched for
#' the optimal threshold when \code{t0} is left unspecified}
#'   \item{\code{t0.plot}}{controls whether the objective function will be plotted.
#' When \code{t0.plot} is ture, both the raw \code{t0.plot} values and the smoothed estimates
#' (using local polynomial regression fitting) are plotted.}
#' }
#'
#' @importFrom stats approx coef ecdf lm model.extract model.matrix optimize printCoefmat
#' @importFrom stats loess.smooth pnorm qnorm runif sd
#' @importFrom survival survfit coxph Surv
#' @importFrom graphics abline legend lines mtext par plot title
#' 
#' @export
#' @references Qian, J., Chiou, S.H., Maye, J.E., Atem, F., Johnson, K.A. and Betensky, R.A. (2018)
#' Threshold regression to accommodate a censored covariate, \emph{Biometrics}, \bold{74}(4):
#' 1261--1270.
#' @references Atem, F., Qian, J., Maye J.E., Johnson, K.A. and Betensky, R.A. (2017),
#' Linear regression with a randomly censored covariate: Application to an Alzheimer's study.
#' \emph{Journal of the Royal Statistical Society: Series C}, \bold{66}(2):313--328.
#'
#' @examples
#' simDat <- function(n) {
#'   X <- rexp(n, 3)
#'   Z <- runif(n, 1, 6)
#'   Y <- 0.5 + 0.5 * X - 0.5 * Z + rnorm(n, 0, .75)
#'   cstime <- rexp(n, .75)
#'   delta <- (X <= cstime) * 1
#'   X <- pmin(X, cstime)
#'   data.frame(Y = Y, X = X, Z = Z, delta = delta)
#' }
#'
#' set.seed(0)
#' dat <- simDat(200)
#'
#' ## Falsely assumes all covariates are free of censoring
#' thlm(Y ~ X + Z, data = dat)
#' 
#' ## Complete cases regression
#' thlm(Y ~ X + Z, cens = delta, data = dat)
#' thlm(Y ~ X + Z, data = dat, subset = delta == 1) ## same results
#' 
#' ## reverse survival regression
#' thlm(Y ~ X + Z, cens = delta, data = dat, method = "reverse")
#' 
#' ## threshold regression without bootstrap
#' thlm(Y ~ X + Z, cens = delta, data = dat, method = "dt")
#' thlm(Y ~ X + Z, cens = delta, data = dat, method = "ct", control =
#' list(t0.interval = c(0.2, 0.6), t0.plot = FALSE))
#' 
#' ## threshold regression with bootstrap
#' thlm(Y ~ X + Z, cens = delta, data = dat, method = "dt", B = 100)
#' thlm(Y ~ X + Z, cens = delta, data = dat, method = "ct", B = 100)
#' 
#' ## display all
#' thlm(Y ~ X + Z, cens = delta, data = dat, method = "all", B = 100)

thlm <- function(formula, data, subset, method = "cc", B = 0,
                 x.upplim = NULL, t0 = NULL, control = thlm.control()) {
    Call <- match.call()
    mnames <- c("", "formula", "data", "subset", "cens")
    cnames <- names(Call)
    cnames <- cnames[match(mnames, cnames, 0)]
    mcall <- Call[cnames]
    mcall[[1]] <- as.name("model.frame")
    obj <- eval(mcall, parent.frame())
    y <- as.numeric(model.extract(obj, "response"))
    nSurv <- sapply(1:ncol(obj), function(x) class(obj[,x]))
    covNames <- all.vars(formula)[-c(1, which(nSurv == "Surv") + 1)]
    if (sum(nSurv == "Surv") > 1) stop("The current version allows at most ONE censored covariate.")
    if (sum(nSurv == "Surv") == 0) {
        out <- lm(formula, data = obj)
        method <- "lm"
    }
    if (sum(nSurv == "Surv") == 1) {
        tmp <- do.call(cbind, sapply(1:ncol(obj), function(x) as.matrix(obj[,x])))
        cens <- tmp$status
        u <- tmp$time
        z <- tmp[, attr(tmp, "dimnames")[[2]] == ""][,-1]
        if (method %in% c("cc", "complete cases")) {
            ## out <- lm(formula, subset = cens == 1)
            out <- cc.reg(y, u, cens, z)
            names(out$a2) <- names(out$a2.sd) <- covNames
            method <- "clm"
        }
        if (method %in% c("reverse", "rev", "reverse survival")) {
            out <- rev.surv.reg(y, u, cens, z)
            method <- "rev"
        }
        if (method %in% c("deletion", "deletion threshold", "dt")) {
            op2 <- par(mar = c(3.5, 3.5, 2.5, 2.5))
            out <- threshold.reg.m1(y, u, cens, z, x.upplim, t0, control)
            names(out$a2) <- names(out$a2.sd) <- covNames
            method <- "dt"
            out$a1.sd <- NA
            par(op2)
        }
        if (method %in% c("complete", "complete threshold", "ct")) {
            op2 <- par(mar = c(3.5, 3.5, 2.5, 2.5))
            out <- threshold.reg.m2(y, u, cens, z, x.upplim, t0, control)
            names(out$a2) <- names(out$a2.sd) <- covNames
            method <- "ct"
            out$a1.sd <- NA
            par(op2)
        }
        if (method == "all") {
            out <- NULL
            out$cc <- cc.reg(y, u, cens, z)
            names(out$cc$a2) <- names(out$cc$a2.sd) <- covNames
            out$rev <- rev.surv.reg(y, u, cens, z)
            op2 <- par(mfrow = c(2, 1), mar = c(3.5, 3.5, 2.5, 2.5))
            out$dt <- threshold.reg.m1(y, u, cens, z, x.upplim, t0, control)
            names(out$dt$a2) <- names(out$dt$a2.sd) <- covNames
            out$dt$a1.sd <- NA
            out$ct <- threshold.reg.m2(y, u, cens, z, x.upplim, t0, control)
            names(out$ct$a2) <- names(out$ct$a2.sd) <- covNames
            out$ct$a1.sd <- NA
            par(op2)
        }
    }
    if (B > 0 & method %in% c("dt", "ct", "all")) {
        a1.bp <- a1.bp2 <- rep(NA, B)
        k <- 1
        t0 <- out$threshold
        n <- length(y)
        while(k <= B) {
            bp.id <- sample(1:n, n, replace = TRUE)
            y.bp <- y[bp.id]
            u.bp <- u[bp.id]
            z.bp <- as.matrix(z[bp.id,])
            delta.bp <- cens[bp.id]
            temp <- diff(sort(u.bp))
            temp2 <- sort(temp[temp>0])
            u.bp <- u.bp + runif(n, 0, min(temp2) / n)
            if((sum(u.bp[delta.bp == 1] <= t0) > 0) & (sum(u.bp > t0) > 0)) {
                if (method == "dt") 
                    a1.bp[k] <- threshold.reg.m1(y.bp, u.bp, delta.bp, z.bp, x.upplim, t0)$a1
                if (method == "ct")
                    a1.bp[k] <- threshold.reg.m2(y.bp, u.bp, delta.bp, z.bp, x.upplim, t0)$a1
                k <- k + 1
            }
            if (method == "all" & (sum(u.bp[delta.bp == 1] <= out$dt$threshold) > 0) & (sum(u.bp > out$dt$threshold) > 0) &
                (sum(u.bp[delta.bp == 1] <= out$ct$threshold) > 0) & (sum(u.bp > out$ct$threshold) > 0)) {
                a1.bp[k] <- threshold.reg.m1(y.bp, u.bp, delta.bp, z.bp, x.upplim, out$dt$threshold)$a1
                a1.bp2[k] <- threshold.reg.m2(y.bp, u.bp, delta.bp, z.bp, x.upplim, out$ct$threshold)$a1
                k <- k + 1
            }
        }
        if (method == "all") {
            out$dt$a1.sd <- sd(a1.bp)
            out$ct$a1.sd <- sd(a1.bp2)
        } else { out$a1.sd <- sd(a1.bp) }
    }
    out$method <- method
    out$Call <- Call
    out$covNames <- names(obj)[-ncol(obj)]
    out$B <- B
    class(out) <- "thlm"
    return(out)
}

cc.reg <- function(y, u, delta, z = NULL) {
    fit <- lm(y ~ u + z, subset = delta == 1)
    alpha1.est <- fit$coef["u"]
    alpha1.sd <- coef(summary(fit))[2,2]
    alpha2.est <- fit$coef[-(1:2)]
    alpha2.sd <- coef(summary(fit))[-(1:2),2]
    power <- (((alpha1.est + alpha1.sd * qnorm(0.975))<0) ||
              (0<(alpha1.est - alpha1.sd * qnorm(0.975))))*1
    ## output results
    list(a1 = alpha1.est, a1.sd = alpha1.sd,
         a2 = alpha2.est, a2.sd = alpha2.sd,
         power = power)
}

rev.surv.reg <- function(y, u, delta, z = NULL) {
    if (is.null(z)) fit <- coxph(Surv(u, delta) ~ y)
    else fit <- coxph(Surv(u, delta) ~ y + z)
    alpha1.est <- fit$coef["y"]
    alpha1.sd <- coef(summary(fit))[1,3]
    alpha1.pval <- coef(summary(fit))[1,5]
    power <- (alpha1.pval < 0.05)*1
    list(a1 = alpha1.est, a1.sd = alpha1.sd, power = power)
}

thlm.control <- function(t0.interval = NULL, t0.plot = TRUE) {
    list(t0.interval = t0.interval, t0.plot = t0.plot)
}

opt.threshold.m1 <- function(y, u, delta, x.upplim = 1.75, t0.interval = NULL, t0.plot = FALSE) {
    n <- length(y)
    if (is.null(t0.interval)) t0.interval <- range(u)
    km.fit <- survfit(Surv(u, delta) ~ 1)
    if (delta[which.max(u)]) {
        km.time <- summary(km.fit)$time
        km.est <- summary(km.fit)$surv
    } else {
        km.time <- c(summary(km.fit)$time, max(u))
        km.est <- c(summary(km.fit)$surv, 0)
    }
    km.time[which.max(km.time)] <- max(km.time, x.upplim)
    obj.m1 <- function(t0){
        ind <- sum(km.time <= t0)
        surv.t0 <- km.est[ind+1] + (km.est[ind] - km.est[ind + 1]) *
            (km.time[ind + 1] - t0) / (km.time[ind + 1] - km.time[ind])
        ## approx(km.time, km.est, t0)$y
        km.time.2 <- c(0, km.time[-length(km.time)])
        km.est.2 <- c(1, km.est[-length(km.est)]) 
        vec.a <- (km.time > t0) * (km.time - pmax(km.time.2, t0)) *
            (pmin(km.est.2, surv.t0) + pmin(km.est, surv.t0)) / 2
        est.cond.x <- sum(vec.a) / surv.t0 + t0
        bias.correct <- est.cond.x - sum(u * delta * (u <= t0)) / sum(delta * (u <= t0))
        n1 <- sum((u <= t0 & delta == 1))
        n2 <- sum(u > t0) 
        return(abs(bias.correct) / (sqrt(1 / n1 + 1 / n2))) ## assuming equal variance    
    }
    t0.opt <- optimize(obj.m1, t0.interval, tol = 1e-5, maximum = TRUE)
    if (t0.plot) {
        t0.vec <- seq(t0.interval[1], t0.interval[2], length = 50)[2:49]
        t0.thr <- unlist(sapply(t0.vec, obj.m1))
        t0.vec <- c(t0.vec, t0.opt$maximum)
        t0.thr <- c(t0.thr, t0.opt$objective)
        t0.thr <- t0.thr[order(t0.vec)]
        t0.vec <- t0.vec[order(t0.vec)]
        plot(t0.vec, t0.thr, "l", lty = 2, lwd = 2, main = "", xlab = "", ylab = "")
        mtext(expression(bold("Threshold estimation with deletion threshold regression")), 3, line = .5, cex = 1.2)
        title(xlab = "Time", ylab = "Objective function", line = 2)
        fit <- loess.smooth(t0.vec, t0.thr, degree = 2)
        lines(fit$x, fit$y, lty = 1, lwd = 2)
        abline(v = t0.opt$maximum, lty = "dotted", lwd = 1.5)
        legend("topright", c("raw value", "smoothed"), lty = 1:2, bty = "n")
    }
    list(t0.opt = t0.opt$maximum, obj.val = t0.opt$objective)
}

threshold.reg.m1 <- function(y, u, delta, z = NULL, x.upplim = 1.75, t0 = NULL, control = thlm.control()) {
    if (is.null(x.upplim)) x.upplim <- max(u)
    if (is.null(t0)) t0 <- opt.threshold.m1(y, u, delta, x.upplim, control$t0.interval, control$t0.plot)$t0.opt
    n <- length(y)
    km.fit <- survfit(Surv(u, delta) ~ 1)
    if (delta[which.max(u)]) {
        km.time <- summary(km.fit)$time
        km.est <- summary(km.fit)$surv
    } else {
        km.time <- c(summary(km.fit)$time, max(u))
        km.est <- c(summary(km.fit)$surv, 0)
    }
    km.time[which.max(km.time)] <- max(km.time, x.upplim)
    ind <- sum(km.time <= t0)
    surv.t0 <- km.est[ind+1] + (km.est[ind] - km.est[ind + 1]) *
        (km.time[ind + 1] - t0) / (km.time[ind + 1] - km.time[ind])
    km.time.2 <- c(0, km.time[-length(km.time)])
    km.est.2 <- c(1, km.est[-length(km.est)]) 
    vec.a <- (km.time > t0) * (km.time - pmax(km.time.2, t0)) *
        (pmin(km.est.2, surv.t0) + pmin(km.est, surv.t0)) / 2
    est.cond.x <- sum(vec.a) / surv.t0 + t0
    delta2 <- 1 * (u > t0)
    fit.lm <- lm(y ~ delta2 + z, subset = ((u <= t0 & delta == 1) | (u > t0)))
    beta1.est <- coef(fit.lm)[2]
    beta1.sd <- coef(summary(fit.lm))[2,2]
    power <- as.numeric(1 * (abs(beta1.est) / beta1.sd > qnorm(.975)))
    beta2.est <- coef(fit.lm)[-(1:2)]
    beta2.sd <- coef(summary(fit.lm))[-(1:2), 2]
    bias.correct <- est.cond.x - sum(u * delta * (u <= t0)) / sum(delta * (u <= t0))
    alpha1.est <- beta1.est / bias.correct
    ## percentage of different portions
    percent.del <- sum((1 - delta) * (u <= t0))/n
    percent.below <- sum(delta * (u <= t0))/n
    percent.above <- sum(u > t0) / n
    list(a1 = alpha1.est, b1 = beta1.est, b1.sd = beta1.sd, power.b1 = power,
         a2 = beta2.est, a2.sd = beta2.sd, cond.x = est.cond.x, pdel = percent.del,
         pbelow = percent.below, pabove = percent.above, bias.cor = bias.correct, threshold = t0)
}

opt.threshold.m2 <- function(y, u, delta, x.upplim = 1.75, t0.interval = NULL, t0.plot = FALSE) {
    n <- length(y)
    if (is.null(t0.interval)) t0.interval <- range(u)
    km.fit <- survfit(Surv(u, delta) ~ 1)
    if (delta[which.max(u)]) {
        km.time <- summary(km.fit)$time
        km.est <- summary(km.fit)$surv
    } else {
        km.time <- c(summary(km.fit)$time, max(u))
        km.est <- c(summary(km.fit)$surv, 0)
    }
    km.time[which.max(km.time)] <- max(km.time, x.upplim)
    obj.m2 <- function(t0) {
        ind <- sum(km.time <= t0)
        surv.t0 <- km.est[ind+1] + (km.est[ind] - km.est[ind + 1]) *
            (km.time[ind + 1] - t0) / (km.time[ind + 1] - km.time[ind])
        km.time.2 <- c(0, km.time[-length(km.time)])
        km.est.2 <- c(1, km.est[-length(km.est)])
        vec.a <- (km.time > t0) * (km.time - pmax(km.time.2, t0)) *
            (pmin(km.est.2, surv.t0) + pmin(km.est, surv.t0)) / 2
        est.cond.x <- sum(vec.a) / surv.t0 + t0 
        vec.b <- (km.time - km.time.2) * (km.est.2 + km.est) / 2
        est.mean.x <- sum(vec.b)
        est.prob.u <- approx(u, ecdf(u)(u), t0)$y
        bias.correct <- (est.cond.x - est.mean.x) / est.prob.u
        n1 <- sum(u <= t0)
        n2 <- sum(u > t0)
        return(abs(bias.correct) / (sqrt(1/n1 + 1/n2)))
    }
    t0.opt <- optimize(obj.m2, t0.interval, tol = 1e-5, maximum = TRUE)
    if (t0.plot) {
        t0.vec <- seq(t0.interval[1], t0.interval[2], length = 50)[2:49]
        t0.thr <- unlist(sapply(t0.vec, obj.m2))
        t0.vec <- c(t0.vec, t0.opt$maximum)
        t0.thr <- c(t0.thr, t0.opt$objective)
        t0.thr <- t0.thr[order(t0.vec)]
        t0.vec <- t0.vec[order(t0.vec)]        
        plot(t0.vec, t0.thr, "l", lty = 2, lwd = 2, main = "", xlab = "", ylab = "")
        mtext(expression(bold("Threshold estimation with complete threshold regression")), 3, line = .5, cex = 1.2)
        title(xlab = "Time", ylab = "Objective function", line = 2)
        fit <- loess.smooth(t0.vec, t0.thr, degree = 2)
        lines(fit$x, fit$y, lty = 1, lwd = 2)
        abline(v = t0.opt$maximum, lty = "dotted", lwd = 1.5)
        legend("topright", c("raw value", "smoothed"), lty = 1:2, bty = "n")
    }
    list(t0.opt = t0.opt$maximum, obj.val = t0.opt$objective)
}


threshold.reg.m2 <- function(y, u, delta, z = NULL, x.upplim = 1.75, t0 = NULL, control = thlm.control()) {
    if (is.null(x.upplim)) x.upplim <- max(u)
    if (is.null(t0)) t0 <- opt.threshold.m2(y, u, delta, x.upplim, control$t0.interval, control$t0.plot)$t0.opt
    n <- length(y)
    km.fit <- survfit(Surv(u, delta) ~ 1)
    if (delta[which.max(u)]) {
        km.time <- summary(km.fit)$time
        km.est <- summary(km.fit)$surv
    } else {
        km.time <- c(summary(km.fit)$time, max(u))
        km.est <- c(summary(km.fit)$surv, 0)
    }
    km.time[which.max(km.time)] <- max(km.time, x.upplim)
    ind <- sum(km.time <= t0)
    surv.t0 <- km.est[ind+1] + (km.est[ind] - km.est[ind + 1]) *
        (km.time[ind + 1] - t0) / (km.time[ind + 1] - km.time[ind])
    km.time.2 <- c(0, km.time[-length(km.time)])
    km.est.2 <- c(1, km.est[-length(km.est)]) 
    vec.a <- (km.time > t0) * (km.time - pmax(km.time.2, t0)) *
        (pmin(km.est.2, surv.t0) + pmin(km.est, surv.t0)) / 2
    est.cond.x <- sum(vec.a) / surv.t0 + t0
    vec.b <- (km.time - km.time.2) * (km.est.2 + km.est) / 2
    est.mean.x <- sum(vec.b)
    est.prob.u <- approx(u, ecdf(u)(u), t0)$y
    delta2 <- 1 * (u > t0)
    fit.lm <- lm(y ~ delta2 + z)
    beta1.est <- coef(fit.lm)[2]
    beta1.sd <- coef(summary(fit.lm))[2,2]
    power <- as.numeric(1 * (abs(beta1.est) / beta1.sd > qnorm(.975)))
    beta2.est <- coef(fit.lm)[-(1:2)]
    beta2.sd <- coef(summary(fit.lm))[-(1:2), 2]
    bias.correct <- (est.cond.x - est.mean.x) / est.prob.u
    alpha1.est <- beta1.est / bias.correct
    percent.below <- sum(u <= t0) / n
    percent.above <- sum(u > t0) / n
    list(a1 = alpha1.est, b1 = beta1.est, b1.sd = beta1.sd, power.b1 = power,
         a2 = beta2.est, a2.sd = beta2.sd, cond.x = est.cond.x, mean.x = est.mean.x,
         pbelow = percent.below, pabove = percent.above, bias.cor = bias.correct, threshold = t0)
}
