#' @export
summary.thlm <- function(object, ...) {
    if (class(object) != "thlm") stop("Must be a thlm object")
    if (object$method == "lm") {
        class(object) <- "lm"
        ## cat(paste("\nLinear regression\n"))
        tab <- coef(summary(object))
        colnames(tab) <- c("Estimate", "StdErr", "z.value", "p.value")
        tab <- round(tab, 7)
    }
    if (object$method == "clm") {
        ## cat(paste("\nComplete-case regression\n"))
        tab <- as.data.frame(cbind(Estimate = c(object$a1, object$a2),
                                   StdErr = c(object$a1.sd, object$a2.sd)))
        tab$z.value <- tab$Estimate / tab$StdErr
        tab$p.value <- with(tab, 2 * pnorm(-abs(Estimate) / StdErr))
        if (any(is.na(object$cc))) tab$cc <- tab$cc[complete.cases(tab$cc),]
        rownames(tab) <- object$names[-1]
        tab <- round(tab, 7)
    }
    if (object$method == "rev") {
        ## cat(paste("\nReverse survival regression\n"))
        tab <- as.data.frame(cbind(Estimate = object$a1, StdErr = object$a1.sd,
                                   z.value = object$a1 / object$a1.sd,
                                   p.value = 2 * pnorm(-abs(object$a1) / object$a1.sd)))
        rownames(tab) <- object$names[1]
        tab <- round(tab, 7)
    }
    if (object$method %in% c("dt", "ct")) {
        if (length(object$a2.sd) == 0) {
            tab <- data.frame(Estimate = c(object$a1, object$b1),
                              StdErr = c(object$a1.sd, object$b1.sd))
            rownames(tab) <- c(paste("a1:", object$names[2], sep = ""),
                               paste("b1:", object$names[2], sep = ""))
        } else {            
            tab <- data.frame(Estimate = c(object$a1, object$a2, object$b1),
                              StdErr = c(object$a1.sd, object$a2.sd, object$b1.sd))
            rownames(tab) <- c(paste("a1:", object$names[2], sep = ""),
                               sapply(3:length(object$names),
                                      function(y) paste("a", y-1, ":", object$names[y], sep = "")),
                               paste("b1:", object$names[2], sep = ""))
        }
        tab$z.value <- tab$Estimate / tab$StdErr
        tab$p.value <- with(tab, 2 * pnorm(-abs(Estimate) / StdErr))
        tab <- round(tab, 7)
    }
    if (object$method == "all") {
        ## complete case
        tab <- NULL
        tab$cc <- as.data.frame(cbind(Estimate = c(object$cc$a1, object$cc$a2),
                                      StdErr = c(object$cc$a1.sd, object$cc$a2.sd)))
        tab$cc$z.value <- tab$cc$Estimate / tab$cc$StdErr
        tab$cc$p.value <- with(tab$cc, 2 * pnorm(-abs(Estimate) / StdErr))
        if (any(is.na(object$cc))) tab$cc <- tab$cc[complete.cases(tab$cc),]
        rownames(tab$cc) <- object$names[-1]
        ## reverse survival regression
        tab$rev <- as.data.frame(cbind(Estimate = object$rev$a1, StdErr = object$rev$a1.sd,
                                       z.value = object$rev$a1 / object$rev$a1.sd,
                                       p.value = 2 * pnorm(-abs(object$rev$a1) / object$rev$a1.sd)))
        rownames(tab$rev) <- object$names[1]
        ## deletion threshold & complete threshold
        if (length(object$dt$a2.sd) == 0) {
            tab$dt <- data.frame(Estimate = c(object$dt$a1, object$dt$b1),
                                 StdErr = c(object$dt$a1.sd, object$dt$b1.sd))
            rownames(tab$dt) <- c(paste("a1:", object$names[2], sep = ""),
                                  paste("b1:", object$names[2], sep = ""))
            tab$ct <- data.frame(Estimate = c(object$ct$a1, object$ct$b1),
                                 StdErr = c(object$ct$a1.sd, object$ct$b1.sd))
            rownames(tab$ct) <- c(paste("a1:", object$names[2], sep = ""),
                                  paste("b1:", object$names[2], sep = ""))
        } else {            
            tab$dt <- data.frame(Estimate = c(object$dt$a1, object$dt$a2, object$dt$b1),
                                 StdErr = c(object$dt$a1.sd, object$dt$a2.sd, object$dt$b1.sd))
            rownames(tab$dt) <- c(paste("a1:", object$names[2], sep = ""),
                                  sapply(3:length(object$names),
                                         function(y) paste("a", y-1, ":", object$names[y], sep = "")),
                                  paste("b1:", object$names[2], sep = ""))
            tab$ct <- data.frame(Estimate = c(object$ct$a1, object$ct$a2, object$ct$b1),
                                 StdErr = c(object$ct$a1.sd, object$ct$a2.sd, object$ct$b1.sd))
            rownames(tab$ct) <- c(paste("a1:", object$names[2], sep = ""),
                                  sapply(3:length(object$names),
                                         function(y) paste("a", y-1, ":", object$names[y], sep = "")),
                                  paste("b1:", object$names[2], sep = ""))
        }
        tab$dt$z.value <- tab$dt$Estimate / tab$dt$StdErr
        tab$dt$p.value <- with(tab$dt, 2 * pnorm(-abs(Estimate) / StdErr))
        object$threshold$dt <- object$dt$threshold
        tab$ct$z.value <- tab$ct$Estimate / tab$ct$StdErr
        tab$ct$p.value <- with(tab$ct, 2 * pnorm(-abs(Estimate) / StdErr))
        object$threshold$ct <- object$ct$threshold
        tab$cc <- round(tab$cc, 7)
        tab$rev <- round(tab$rev, 7)
        tab$dt <- round(tab$dt, 7)
        tab$ct <- round(tab$ct, 7)
    }
    out <- list(Call = object$Call, method = object$method, tab = tab, threshold = object$threshold)
    class(out) <- "summary.thlm"
    out
}

#' @export
print.summary.thlm <- function(x, ...) {
    cat("Call:\n")
    print(x$Call)
    if (x$method == "lm") cat(paste("\n Linear regression\n"))
    if (x$method == "clm") cat(paste("\n Complete-case regression\n"))
    if (x$method == "rev") cat(paste("\n Reverse survival regression\n"))
    if (x$method == "dt") cat(paste("\n Deletion threshold regression\n"))
    if (x$method == "ct") cat(paste("\n Complete threshold regression\n"))
    if (x$method %in% c("dt", "ct"))
        cat(paste(" Optimal threshold is", round(x$threshold, 3), "\n"))
    if (x$method != "all") {
        cat("\n Coefficients:\n")
        printCoefmat(x$tab, P.values = TRUE, has.Pvalue = TRUE)
        cat("\n")
        if (x$method %in% c("dt", "ct")) {
            cat("Under the assumption that X is independent of Z given X*. See ?thlm for more detail.")
            cat("\n")
        }
    }
    if (x$method == "all") {
        cat(paste("\nComplete-case regression\n"))
        cat("\n Coefficients:\n")
        printCoefmat(x$tab$cc, P.values = TRUE, has.Pvalue = TRUE)
        cat(paste("\nReverse survival regression\n"))
        cat("\n Coefficients:\n")
        printCoefmat(x$tab$rev, P.values = TRUE, has.Pvalue = TRUE)
        cat(paste("\nDeletion threshold regression\n"))
        cat(paste(" Optimal threshold is", round(x$threshold$dt, 3)))
        cat("\n Coefficients:\n")
        printCoefmat(x$tab$dt, P.values = TRUE, has.Pvalue = TRUE)
        cat(paste("\nComplete threshold regression\n"))
        cat(paste(" Optimal threshold is", round(x$threshold$ct, 3)))
        cat("\n Coefficients:\n")
        printCoefmat(x$tab$ct, P.values = TRUE, has.Pvalue = TRUE)
        cat("\n")
        cat("Under the assumption that X is independent of Z given X*. See ?thlm for more detail.")
        cat("\n")
    }
}

#' @export
print.thlm <- function(x, ...) {
    if (class(x) != "thlm") stop("Must be a thlm object")
    cat("\n Call: ")
    print(x$Call)
    if (x$method %in% c("clm", "rev")) {
        cat("\n Hypothesis test of association\n")
        cat(paste("H0: a1 = 0, p-value =",
                  sprintf("%.4f", round(2 * pnorm(-abs(x$a1) / x$a1.sd), 4)), "\n"))
    }
    if (x$method %in% c("dt", "ct")) {
        if (x$B == 0) {
            cat("\n Hypothesis test of association\n")
            cat(paste(" H0: b1 = 0, p-value =", sprintf("%.4f", round(2 * pnorm(-abs(x$b1) / x$b1.sd), 4)), "\n"))
        } else {
            cat("\n Hypothesis test of association\n")
            cat(paste(" H0: b1 = 0, p-value =",
                      sprintf("%.4f", round(2 * pnorm(-abs(x$b1) / x$b1.sd), 4))))
            cat(paste("\n H0: a1 = 0, p-value =",
                      sprintf("%.4f", round(2 * pnorm(-abs(x$a1) / x$a1.sd), 4)), "\n"))
        }
    }
    if (x$method == "lm") {
        cat("\n Hypothesis test of association\n")
        class(x) <- "lm"
        cat(paste(" H0: a1 = 0, p-value =",
                  sprintf("%.4f", round(coef(summary(x))[2,4], 4)), "\n"))
    }
    if (x$method == "all") {
        cat("\n Hypothesis test of association\n")
        cat("\n Complete-cases\n")
        cat(paste(" H0: a1 = 0, p-value =",sprintf("%.4f", round(2 * pnorm(-abs(x$cc$a1) / x$cc$a1.sd), 4)), "\n"))
        cat("\n Reverse survival\n")
        cat(paste(" H0: a1 = 0, p-value =",sprintf("%.4f", round(2 * pnorm(-abs(x$rev$a1) / x$rev$a1.sd), 4)), "\n"))        
        if (x$B == 0) {
            cat("\n Deletion threshold\n")
            cat(paste(" H0: b1 = 0, p-value =",sprintf("%.4f", round(2 * pnorm(-abs(x$dt$b1) / x$dt$b1.sd), 4)), "\n"))
            cat("\n Complete threshold\n")
            cat(paste(" H0: b1 = 0, p-value =",sprintf("%.4f", round(2 * pnorm(-abs(x$ct$b1) / x$ct$b1.sd), 4)), "\n"))
        } else {
            cat("\n Deletion threshold\n")
            cat(paste(" H0: b1 = 0: p-value =", sprintf("%.4f", round(2 * pnorm(-abs(x$dt$b1) / x$dt$b1.sd),4))))
            cat(paste("\n H0: a1 = 0: p-value =", sprintf("%.4f", round(2 * pnorm(-abs(x$dt$a1) / x$dt$a1.sd),4))))
            cat("\n\n Complete threshold\n")
            cat(paste(" H0: b1 = 0: p-value =", sprintf("%.4f", round(2 * pnorm(-abs(x$ct$b1) / x$ct$b1.sd),4))))
            cat(paste("\n H0: a1 = 0: p-value =", sprintf("%.4f", round(2 * pnorm(-abs(x$ct$a1) / x$ct$a1.sd),4))))
        }
    }
    cat("\n")
}
