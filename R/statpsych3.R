# ======================== Confidence Intervals ==============================
#  ci.prop1 ================================================================== 
#' Confidence intervals for a single proportion
#'
#'
#' @description
#' Computes adjusted Wald and Wilson confidence intervals for a single
#' population proportion. The Wilson confidence interval uses a continuity
#' correction.
#'
#'
#' @param   alpha   alpha level for 1-alpha confidence
#' @param   f       number of participants who have the attribute
#' @param   n       sample size
#'
#'
#' @return
#' Returns a 2-row matrix. The columns of row 1 are:
#' * Estimate - adjusted estimate of proportion
#' * SE - adjusted standard error
#' * LL - lower limit of the adjusted Wald confidence interval
#' * UL - upper limit of the adjusted Wald confidence interval
#'
#'
#' The columns of row 2 are:
#' * Estimate - ML estimate of proportion
#' * SE - standard error
#' * LL - lower limit of the Wilson confidence interval
#' * UL - upper limit of the Wilson confidence interval
#'
#'
#' @references
#' \insertRef{Agresti1998}{statpsych}
#'
#'
#' @examples
#' ci.prop1(.05, 12, 100)
#'
#' # Should return:
#' #                  Estimate         SE         LL        UL
#' # Adjusted Wald   0.1346154 0.03346842 0.06901848 0.2002123
#' # Wilson with cc  0.1200000 0.03249615 0.06625153 0.2039772
#'
#'
#' @importFrom stats qnorm
#' @export
ci.prop1 <- function(alpha, f, n) {
 if (f > n) {stop("f cannot be greater than n")}
 z <- qnorm(1 - alpha/2)
 p.mle <- f/n
 se.mle <- sqrt(p.mle*(1 - p.mle)/n)
 b1 <- 2*n*p.mle + z^2
 b2 <- 2*(n + z^2)
 LL.wil <- (b1 - 1 - z*sqrt(z^2 - 2 - 1/n + 4*p.mle*(n*(1 - p.mle) + 1)))/b2
 UL.wil <- (b1 + 1 + z*sqrt(z^2 + 2 - 1/n + 4*p.mle*(n*(1 - p.mle) - 1)))/b2
 if (p.mle == 0) {LL.wil = 0}
 if (p.mle == 1) {UL.wil = 1}
 p.adj <- (f + 2)/(n + 4)
 se.adj <- sqrt(p.adj*(1 - p.adj)/(n + 4))
 LL.adj <- p.adj - z*se.adj
 UL.adj <- p.adj + z*se.adj
 if (LL.adj < 0) {LL.adj = 0}
 if (UL.adj > 1) {UL.adj = 1}
 out1 <- t(c(p.adj, se.adj, LL.adj, UL.adj))
 out2 <- t(c(p.mle, se.mle, LL.wil, UL.wil))
 out <- rbind(out1, out2)
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 rownames(out) <- c("Adjusted Wald", "Wilson with cc")
 return(out)
}


#  ci.prop.fpc =============================================================== 
#' Confidence interval for a single proportion with a finite population 
#' correction
#'
#'
#' @description
#' Computes an adjusted Wald interval for a single population proportion with a 
#' finite population correction (fpc). This confidence interval is useful 
#' when the sample size is not a small fraction of the population size.
#'
#'
#' @param   alpha   alpha level for 1-alpha confidence
#' @param   f       number of participants who have the attribute
#' @param   n       sample size
#' @param   N       population size
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - adjusted estimate of proportion
#' * SE - adjusted standard error with fpc
#' * LL - lower limit of the confidence interval with fpc
#' * UL - upper limit of the confidence interval with fpc
#'
#'
#' @examples
#' ci.prop.fpc(.05, 12, 100, 400)
#'
#' # Should return:
#' #   Estimate         SE         LL        UL
#' #  0.1346154  0.0290208 0.07773565 0.1914951
#'
#'
#' @importFrom stats qnorm
#' @export
ci.prop.fpc <- function(alpha, f, n, N) {
 if (f > n) {stop("f cannot be greater than n")}
 if (f > N) {stop("f cannot be greater than N")}
 if (n > N) {stop("n cannot be greater than N")}
 z <- qnorm(1 - alpha/2)
 p.adj <- (f + 2)/(n + 4)
 se.adj <- sqrt(p.adj*(1 - p.adj)/(n + 4))*sqrt((N - n)/(N - 1))
 LL.adj <- p.adj - z*se.adj
 UL.adj <- p.adj + z*se.adj
 if (LL.adj < 0) {LL.adj = 0}
 if (UL.adj > 1) {UL.adj = 1}
 out <- t(c(p.adj, se.adj, LL.adj, UL.adj))
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 rownames(out) <- ""
 return(out)
}


#  ci.pairs.mult ============================================================
#' Confidence intervals for pairwise proportion differences of a
#' multinomial variable in a single sample
#'
#'
#' @description
#' Computes adjusted Wald confidence intervals for pairwise proportion
#' differences of a multinomial variable in a single sample. These adjusted
#' Wald confidence intervals use the same method that is used to compare the
#' two proportions in a paired-samples design.
#'
#'
#' @param   alpha   alpha level for 1-alpha confidence
#' @param   f       vector of multinomial frequency counts 
#'
#'
#' @return
#' Returns a matrix with the number of rows equal to the number
#' of pairwise comparisons. The columns are:
#' * Estimate - adjusted estimate of proportion difference
#' * SE - adjusted standard error
#' * LL - lower limit of the adjusted Wald confidence interval
#' * UL - upper limit of the adjusted Wald confidence interval
#'
#'
#' @references
#' \insertRef{Bonett2012}{statpsych}
#'
#'
#' @examples
#' f <- c(125, 82, 92)
#' ci.pairs.mult(.05, f)
#'
#' # Should return:
#' #        Estimate         SE          LL         UL
#' # 1 2  0.14285714 0.04731825  0.05011508 0.23559920
#' # 1 3  0.10963455 0.04875715  0.01407230 0.20519680
#' # 2 3 -0.03322259 0.04403313 -0.11952594 0.05308076
#'
#'
#' @importFrom stats qnorm
#' @export
ci.pairs.mult <-function(alpha, f) {
 zcrit <- qnorm(1 - alpha/2)
 a <- length(f)
 n <- sum(f)
 p.ml <- f/n
 diff.ml <- outer(p.ml, p.ml, '-')
 diff.ml <- (-1)*diff.ml[lower.tri(diff.ml)]
 p <- (f + 1)/(n + 2)
 diff <- outer(p, p, '-')
 diff <- (-1)*diff[lower.tri(diff)]
 v1 <- outer(p, p, "+")
 v2 <- (outer(p, p, "-"))^2
 v <- (v1 - v2)/(n + 2)
 SE <- sqrt(v[lower.tri(v)])
 LL <- diff - zcrit*SE
 UL <- diff + zcrit*SE
 Estimate <- diff
 pair <- t(combn(seq(1:a), 2))
 out <- cbind(pair, Estimate, SE, LL, UL)
 rownames(out) <- rep("", a*(a - 1)/2)
 return(out)
}


#  ci.prop.inv =============================================================== 
#' Confidence interval for a single proportion using inverse sampling
#'
#'
#' @description
#' Computes an exact confidence interval for a single population proportion 
#' when inverse sampling has been used. An approximate standard error is 
#' recovered from the confidence interval. With inverse sampling, the number 
#' of participants who have the attribute (f) is predetermined and sampling 
#' continues until f attains its prespecified value. With inverse sampling, 
#' the sample size (n) will not be known in advance.
#'
#'
#' @param   alpha   alpha level for 1-alpha confidence
#' @param   f       number of participants who have the attribute (fixed)
#' @param   n       sample size (random)
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimate of proportion
#' * SE - recovered standard error
#' * LL - lower limit of confidence interval
#' * UL - upper limit of confidence interval
#'
#'
#' @references
#' \insertRef{Zou2010}{statpsych}
#'
#'
#' @examples
#' ci.prop.inv(.05, 5, 67)
#'
#' # Should return:
#' #   Estimate         SE         LL        UL
#' # 0.07462687 0.03145284 0.02467471 0.1479676
#'
#'
#' @importFrom stats qnorm
#' @importFrom stats qf
#' @export
ci.prop.inv <- function(alpha, f, n) {
 if (f > n) {stop("f cannot be greater than n")}
 z <- qnorm(1 - alpha/2)
 y <- n - f
 est <- f/n
 df1 <- 2*(y + 1)
 df2 <- 2*f
 df3 <- 2*y
 fcritL <- qf(1 - alpha/2, df1, df2)
 ll <- 1/(1 + fcritL*(y + 1)/f)
 if (y == 0) {
   ul <- 1
 } else {
   fcritU <- qf(alpha/2, df3, df2)
   ul <- 1/(1 + fcritU*(y/f))
 }
 se <- (ul - ll)/(2*z)
 out <- t(c(est, se, ll, ul))
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 rownames(out) <- ""
 return(out)
}


#  ci.prop2 ==================================================================
#' Confidence interval for a 2-group proportion difference
#'
#' @description
#' Computes an adjusted Wald confidence interval for a proportion difference
#' in a 2-group design.
#'
#'
#' @param   alpha   alpha level for 1-alpha confidence
#' @param   f1      number of participants in group 1 who have the attribute
#' @param   f2      number of participants in group 2 who have the attribute
#' @param   n1      sample size for group 1
#' @param   n2      sample size for group 2
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - adjusted estimate of proportion difference
#' * SE - adjusted standard error 
#' * LL - lower limit of the adjusted Wald confidence interval
#' * UL - upper limit of the adjusted Wald confidence interval
#'
#'
#' @references
#' \insertRef{Agresti2000}{statpsych}
#'
#'
#' @examples
#' ci.prop2(.05, 35, 21, 150, 150)
#'
#' # Should return:
#' #   Estimate         SE          LL        UL
#' # 0.09210526 0.04476077 0.004375769 0.1798348
#'
#'
#' @importFrom stats qnorm
#' @export
ci.prop2 <- function(alpha, f1, f2, n1, n2) {
 if (f1 > n1) {stop("f cannot be greater than n")}
 if (f2 > n2) {stop("f cannot be greater than n")}
 z <- qnorm(1 - alpha/2)
 p1 <- (f1 + 1)/(n1 + 2)
 p2 <- (f2 + 1)/(n2 + 2)
 se <- sqrt(p1*(1 - p1)/(n1 + 2) + p2*(1 - p2)/(n2 + 2))
 ll <- p1 - p2 - z*se
 ul <- p1 - p2 + z*se
 out <- t(c(p1-p2, se, ll, ul))
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 rownames(out) <- ""
 return(out)
}


# ci.prop2.inv ================================================================
#' Confidence interval for a 2-group proportion difference using inverse 
#' sampling
#'
#'
#' @description
#' Computes an approximate confidence interval for a population proportion 
#' difference when inverse sampling has been used. An approximate standard  
#' error is recovered from the confidence interval. With inverse sampling, the  
#' number of participants who have the attribute within group 1 (f1) and group 2
#' (f2) are predetermined, and sampling continues within each group until f1 
#' and f2 attain their prespecified values. With inverse sampling, the sample 
#' sizes (n1 and n2) will not be known in advance.
#'
#'
#' @param   alpha   alpha level for 1-alpha confidence
#' @param   f1      number of participants in group 1 who have the attribute (fixed)
#' @param   f2      number of participants in group 2 who have the attribute (fixed)
#' @param   n1      sample size for group 1 (random)
#' @param   n2      sample size for group 2 (random)
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimate of proportion difference
#' * SE - recovered standard error 
#' * LL - lower limit of confidence interval
#' * UL - upper limit of confidence interval
#'
#'
#' @references
#' \insertRef{Zou2010}{statpsych}
#'
#'
#' @examples
#' ci.prop2.inv(.05, 10, 10, 48, 213)
#'
#' # Should return:
#' #  Estimate         SE         LL        UL
#' #  0.161385 0.05997618 0.05288277 0.2879851
#'
#'
#' @importFrom stats qnorm
#' @importFrom stats qf
#' @export
ci.prop2.inv <- function(alpha, f1, f2, n1, n2) {
 if (f1 > n1) {stop("f cannot be greater than n")}
 if (f2 > n2) {stop("f cannot be greater than n")}
 z <- qnorm(1 - alpha/2)
 y1 <- n1 - f1
 est1 <- f1/n1
 df1 <- 2*(y1 + 1)
 df2 <- 2*f1
 df3 <- 2*y1
 fcritL <- qf(1 - alpha/2, df1, df2)
 ll1 <- 1/(1 + fcritL*(y1 + 1)/f1)
 if (y1 == 0) {
   ul1 <- 1
 } else {
   fcritU <- qf(alpha/2, df3, df2)
   ul1 <- 1/(1 + fcritU*(y1/f1))
 }
 y2 <- n2 - f2
 est2 <- f2/n2
 df1 <- 2*(y2 + 1)
 df2 <- 2*f2
 df3 <- 2*y2
 fcritL <- qf(1 - alpha/2, df1, df2)
 ll2 <- 1/(1 + fcritL*(y2 + 1)/f2)
 if (y2 == 0) {
   ul2 <- 1
 } else {
   fcritU <- qf(alpha/2, df3, df2)
   ul2 <- 1/(1 + fcritU*(y2/f2))
 }
 diff <- est1 - est2
 ll <- diff - sqrt((est1 - ll1)^2 + (ul2 - est2)^2)
 ul <- diff + sqrt((ul1 - est1)^2 + (est2 - ll2)^2)
 se <- (ul - ll)/(2*z)
 out <- t(c(diff, se, ll, ul))
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 rownames(out) <- ""
 return(out)
}


#  ci.ratio.prop2 ============================================================
#' Confidence interval for a 2-group proportion ratio
#'
#'
#' @description
#' Computes an adjusted Wald confidence interval for a proportion ratio in a
#' 2-group design.
#'
#'
#' @param   alpha   alpha level for 1-alpha confidence
#' @param   f1      number of participants in group 1 who have the attribute
#' @param   f2      number of participants in group 2 who have the attribute
#' @param   n1      sample size for group 1
#' @param   n2      sample size for group 2
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - adjusted estimate of proportion ratio
#' * LL - lower limit of the adjusted Wald confidence interval
#' * UL - upper limit of the adjusted Wald confidence interval
#'
#'
#' @references
#' \insertRef{Price2008}{statpsych}
#'
#'
#' @examples
#' ci.ratio.prop2(.05, 35, 21, 150, 150)
#'
#' # Should return:
#' # Estimate       LL       UL
#' # 1.666667 1.017253 2.705025
#'
#'
#' @importFrom stats qnorm
#' @export
ci.ratio.prop2 <- function(alpha, f1, f2, n1, n2) {
 if (f1 > n1) {stop("f cannot be greater than n")}
 if (f2 > n2) {stop("f cannot be greater than n")}
 z <- qnorm(1 - alpha/2)
 p1 <- (f1 + 1/4)/(n1 + 7/4)
 p2 <- (f2 + 1/4)/(n2 + 7/4)
 v1 <- 1/(f1 + 1/4 + (f1 + 1/4)^2/(n1 - f1 + 3/2))
 v2 <- 1/(f2 + 1/4 + (f2 + 1/4)^2/(n2 - f2 + 3/2))
 se <- sqrt(v1 + v2)
 ll <- exp(log(p1/p2) - z*se)
 ul <- exp(log(p1/p2) + z*se)
 out <- t(c(p1/p2, ll, ul))
 colnames(out) <- c("Estimate", "LL", "UL")
 rownames(out) <- ""
 return(out)
}


#  ci.lc.prop.bs =============================================================
#' Confidence interval for a linear contrast of proportions in a between-
#' subjects design
#'
#'
#' @description
#' Computes an adjusted Wald confidence interval for a linear contrast of 
#' proportions in a betwen-subjects design.
#'
#'
#' @param   alpha   alpha level for 1-alpha confidence
#' @param   f       vector of frequency counts of participants who have the attribute
#' @param   n       vector of sample sizes
#' @param   v       vector of between-subjects contrast coefficients
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - adjusted estimate of proportion linear contrast
#' * SE - adjusted standard error
#' * z - z test statistic
#' * p - two-sided p-value
#' * LL - lower limit of the adjusted Wald confidence interval
#' * UL - upper limit of the adjusted Wald confidence interval
#'
#'
#' @references
#' \insertRef{Price2004}{statpsych}
#'
#'
#' @examples
#' f <- c(26, 24, 38)
#' n <- c(60, 60, 60)
#' v <- c(-.5, -.5, 1)
#' ci.lc.prop.bs(.05, f, n, v)
#'
#' # Should return:
#' #  Estimate         SE        z           p         LL        UL
#' # 0.2119565 0.07602892 2.787841 0.005306059 0.06294259 0.3609705
#'
#'
#' @importFrom stats qnorm
#' @importFrom stats pnorm
#' @export
ci.lc.prop.bs <- function(alpha, f, n, v) {
 s <- sum(as.integer(as.logical(n < f)))
 if (s > 0) {stop("f cannot be greater than n")}
 z <- qnorm(1 - alpha/2)
 m <- length(v) - length(which(v==0))
 p <- (f + 2/m)/(n + 4/m)
 est <- t(v)%*%p
 se <- sqrt(t(v)%*%diag(p*(1 - p))%*%solve(diag(n + 4/m))%*%v)
 zval <- est/se
 pval <- 2*(1 - pnorm(abs(zval)))
 ll <- est - z*se
 ul <- est + z*se
 out <- t(c(est, se, zval, pval, ll, ul))
 colnames(out) <- c("Estimate", "SE", "z", "p", "LL", "UL")
 rownames(out) <- ""
 return(out)
}


#  ci.pairs.prop.bs ========================================================== 
#' Bonferroni confidence intervals for all pairwise proportion differences
#' in a between-subjects design
#'
#'
#' @description
#' Computes adjusted Wald confidence intervals for all pairwise differences of 
#' proportions in a between-subjects design with a Bonferroni adjusted alpha
#' level.
#'
#'
#' @param   alpha   alpha level for simultaneous 1-alpha confidence
#' @param   f       vector of frequency counts of participants who have the attribute
#' @param   n       vector of sample sizes
#'
#'
#' @return
#' Returns a matrix with the number of rows equal to the number
#' of pairwise comparisons. The columns are:
#' * Estimate - adjusted estimate of proportion difference
#' * SE - adjusted standard error
#' * z - z test statistic
#' * p - two-sided p-value
#' * LL - lower limit of the adjusted Wald confidence interval
#' * UL - upper limit of the adjusted Wald confidence interval
#'
#'
#' @references
#' \insertRef{Agresti2000}{statpsych}
#'
#'
#' @examples
#' f <- c(111, 161, 132)
#' n <- c(200, 200, 200)
#' ci.pairs.prop.bs(.05, f, n)
#'
#' # Should return:
#' #        Estimate         SE         z            p          LL          UL
#' # 1 2  -0.2475248 0.04482323 -5.522243 3.346989e-08 -0.35483065 -0.14021885
#' # 1 3  -0.1039604 0.04833562 -2.150803 3.149174e-02 -0.21967489  0.01175409
#' # 2 3   0.1435644 0.04358401  3.293968 9.878366e-04  0.03922511  0.24790360
#'
#'
#' @importFrom stats qnorm
#' @importFrom stats pnorm
#' @export
ci.pairs.prop.bs <-function(alpha, f, n) {
 s <- sum(as.integer(as.logical(n < f)))
 if (s > 0) {stop("f cannot be greater than n")}
 a <- length(f)
 adjalpha <- 2*alpha/(a*(a - 1))
 zcrit <- qnorm(1 - adjalpha/2)
 p.adj <- (f + 1)/(n + 2)
 v <- p.adj*(1 - p.adj)/(n + 2)
 p.adj <- outer(p.adj, p.adj, '-')
 Estimate <- p.adj[upper.tri(p.adj)]
 v <- outer(v, v, "+")
 SE <- sqrt(v[upper.tri(v)])
 z <- Estimate/SE
 p <- 2*(1 - pnorm(abs(z)))
 LL <- Estimate - zcrit*SE
 UL <- Estimate + zcrit*SE
 pair <- t(combn(seq(1:a), 2))
 out <- cbind(pair, Estimate, SE, z, p, LL, UL)
 rownames(out) <- rep("", a*(a - 1)/2)
 return(out)
}


#  ci.slope.prop.bs ========================================================== 
#' Confidence interval for a slope of a proportion in a single-factor 
#' experimental design with a quantitative between-subjects factor
#'
#'
#' @description
#' Computes a test statistic and an adjusted Wald confidence interval for the 
#' slope of proportions in a single-factor experimental design with a 
#' quantitative between-subjects factor. 
#'
#'
#' @param   alpha   alpha level for 1-alpha confidence
#' @param   f       vector of frequency counts of participants who have the attribute
#' @param   n       vector of sample sizes
#' @param   x       vector of quantitative factor values
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - adjusted slope estimate
#' * SE - adjusted standard error
#' * z - z test statistic
#' * p - two-sided p-value
#' * LL - lower limit of the adjusted Wald confidence interval
#' * UL - upper limit of the adjusted Wald confidence interval
#'
#'
#' @references
#' \insertRef{Price2004}{statpsych}
#'
#'
#' @examples
#' f <- c(14, 27, 38)
#' n <- c(100, 100, 100)
#' x <- c(10, 20, 40)
#' ci.slope.prop.bs(.05, f, n, x)
#'
#' # Should return:
#' #    Estimate          SE        z           p          LL         UL
#' # 0.007542293 0.002016793 3.739746 0.000184206 0.003589452 0.01149513
#'
#'
#' @importFrom stats qnorm
#' @importFrom stats pnorm
#' @export
ci.slope.prop.bs <- function(alpha, f, n, x) {
 s <- sum(as.integer(as.logical(n < f)))
 if (s > 0) {stop("f cannot be greater than n")}
 z <- qnorm(1 - alpha/2)
 xmean <- mean(x)
 ssx <- sum((x - xmean)^2)
 c <- (x - xmean)/ssx
 m <- length(c) - length(which(c==0))
 p <- (f + 2/m)/(n + 4/m)
 slope <- t(c)%*%p
 se <- sqrt(t(c)%*%diag(p*(1 - p))%*%solve(diag(n + 4/m))%*%c)
 t <- slope/se
 pval <- 2*(1 - pnorm(abs(t)))
 ll <- slope - z*se
 ul <- slope + z*se
 out <- t(c(slope, se, t, pval, ll, ul))
 colnames(out) <- c("Estimate", "SE", "z", "p", "LL", "UL")
 rownames(out) <- ""
 return(out)
}


#  ci.prop.ps ================================================================
#' Confidence interval for a paired-samples proportion difference
#'
#'
#' @description
#' Computes an adjusted Wald confidence interval for a difference of 
#' proportions in a paired-samples design. This function requires the 
#' frequency counts from a 2 x 2 contingency table for two repeated 
#' dichtomous measurements.
#'
#'
#' @param   alpha  alpha level for 1-alpha confidence
#' @param   f00    number of participants with y = 0 and x = 0
#' @param   f01    number of participants with y = 0 and x = 1
#' @param   f10    number of participants with y = 1 and x = 0
#' @param   f11    number of participants with y = 1 and x = 1
#'
#'
#' @references
#' \insertRef{Bonett2012}{statpsych}
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - adjusted estimate of proportion difference
#' * SE - adjusted standard error
#' * LL - lower limit of the adjusted Wald confidence interval
#' * UL - upper limit of the adjusted Wald confidence interval
#'
#'
#' @examples
#' ci.prop.ps(.05, 12, 4, 26, 6)
#'
#' # Should return:
#' # Estimate         SE        LL        UL
#' #     0.44 0.09448809 0.2548067 0.6251933
#'
#'
#' @importFrom stats qnorm
#' @export
ci.prop.ps <- function(alpha, f00, f01, f10, f11) {
 z <- qnorm(1 - alpha/2)
 n <- f00 + f01 + f10 + f11
 p01 <- (f01 + 1)/(n + 2)
 p10 <- (f10 + 1)/(n + 2)
 diff <- p10 - p01
 se <- sqrt(((p01 + p10) - (p01 - p10)^2)/(n + 2))
 ll <- p10 - p01 - z*se
 ul <- p10 - p01 + z*se
 out <- t(c(diff, se, ll, ul))
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 rownames(out) <- ""
 return(out)
}


#  ci.ratio.prop.ps ==========================================================
#' Confidence interval for a paired-samples proportion ratio
#'
#'
#' @description
#' Computes a confidence interval for a ratio of proportions in a 
#' paired-samples design. This function requires the frequency counts from
#' a 2 x 2 contingency table for two repeated dichotomous measurements.
#'
#'
#' @param   alpha  alpha level for 1-alpha confidence
#' @param   f00    number of participants with y = 0 and x = 0
#' @param   f01    number of participants with y = 0 and x = 1
#' @param   f10    number of participants with y = 1 and x = 0
#' @param   f11    number of participants with y = 1 and x = 1
#'
#'
#' @references
#' \insertRef{Bonett2006a}{statpsych}
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimate of proportion ratio
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @examples
#' ci.ratio.prop.ps(.05, 12, 4, 26, 6)
#'
#' # Should return:
#' # Estimate       LL       UL
#' #      3.2 1.766544 5.796628
#'
#'
#' @importFrom stats qnorm
#' @export
ci.ratio.prop.ps <- function(alpha, f00, f01, f10, f11) {
 z <- qnorm(1 - alpha/2)
 f1 <- f11 + f10
 f2 <- f11 + f01
 n0 <- f11 + f01 + f10
 p1 <- f1/n0
 p2 <- f2/n0
 ratio <- f1/f2
 p1a <- (f1 + 1)/(n0 + 2)
 p2a <- (f2 + 1)/(n0 + 2)
 se.lnp1 <- sqrt((1 - p1a)/((n0 + 2)*p1a)) 
 se.lnp2 <- sqrt((1 - p2a)/((n0 + 2)*p2a))
 se.diff <- sqrt((f10 + f01 + 2)/((f1 + 1)*(f2 + 1)))
 k <- se.diff/(se.lnp1 + se.lnp2)
 z0 <- k*z
 b = 2*(n0 + z0^2)
 LL1 <- (2*f1 + z0^2 - z0*sqrt(z0^2 + 4*f1*(1 - p1)))/b
 UL1 <- (2*f1 + z0^2 + z0*sqrt(z0^2 + 4*f1*(1 - p1)))/b
 LL2 <- (2*f2 + z0^2 - z0*sqrt(z0^2 + 4*f2*(1 - p2)))/b
 UL2 <- (2*f2 + z0^2 + z0*sqrt(z0^2 + 4*f2*(1 - p2)))/b
 ll <- exp(log(LL1) - log(UL2))
 ul <- exp(log(UL1) - log(LL2))
 out <- t(c(ratio, ll, ul))
 colnames(out) <- c("Estimate", "LL", "UL")
 rownames(out) <- ""
 return(out)
}


# ci.condslope.log ===========================================================
#' Confidence intervals for conditional (simple) slopes in a logistic model
#'
#'
#' @description
#' Computes confidence intervals and test statistics for population 
#' conditional slopes (simple slopes) in a logistic model that
#' includes a predictor variable (x1), a moderator variable (x2),
#' and a product predictor variable (x1*x2). Conditional slopes are 
#' computed at low and high values of the moderator variable. 
#'
#'
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  b1     estimated slope coefficient for predictor variable
#' @param  b2     estimated slope coefficient for product variable
#' @param  se1    standard error for predictor coefficient
#' @param  se2    standard error for product coefficient
#' @param  cov    estimated covariance between predictor and product coefficients
#' @param  lo     low value of moderator variable 
#' @param  hi     high value of moderator variable 
#'
#'
#' @return 
#' Returns a 2-row matrix. The columns are:
#' * Estimate - estimated conditional slope
#' * exp(Estimate) - estimated exponentiated conditional slope
#' * z - z test statistic
#' * p - two-sided p-value
#' * LL - lower limit of the exponentiated confidence interval
#' * UL - upper limit of the exponentiated confidence interval
#' 
#' 
#' @examples
#' ci.condslope.log(.05, .132, .154, .031, .021, .015, 5.2, 10.6)
#'
#' # Should return:
#' #                   Estimate exp(Estimate)        z           p 
#' # At low moderator    0.9328      2.541616 2.269824 0.023218266 
#' # At high moderator   1.7644      5.838068 2.906507 0.003654887 
#' #                          LL        UL
#' # At low moderator   1.135802  5.687444
#' # At high moderator  1.776421 19.186357
#'
#'
#' @importFrom stats qnorm
#' @importFrom stats pnorm
#' @export
ci.condslope.log <- function(alpha, b1, b2, se1, se2, cov, lo, hi) {
 z <- qnorm(1 - alpha/2)
 slope.lo <- b1 + b2*lo
 slope.hi <- b1 + b2*hi
 exp.slope.lo <- exp(slope.lo)
 exp.slope.hi <- exp(slope.hi)
 se.lo <- sqrt(se1^2 + se2^2*lo^2 + 2*lo*cov)
 se.hi <- sqrt(se1^2 + se2^2*hi^2 + 2*hi*cov)
 z.lo <- slope.lo/se.lo
 z.hi <- slope.hi/se.hi
 p.lo <- 2*(1 - pnorm(abs(z.lo)))
 p.hi <- 2*(1 - pnorm(abs(z.hi)))
 LL.lo <- exp(slope.lo - z*se.lo)
 UL.lo <- exp(slope.lo + z*se.lo)
 LL.hi <- exp(slope.hi - z*se.hi)
 UL.hi <- exp(slope.hi + z*se.hi)
 out1 <- t(c(slope.lo, exp.slope.lo, z.lo, p.lo, LL.lo, UL.lo))
 out2 <- t(c(slope.hi, exp.slope.hi, z.hi, p.hi, LL.hi, UL.hi))
 out <- rbind(out1, out2)
 colnames(out) <- c("Estimate", "exp(Estimate)", "z", "p", "LL", "UL")
 rownames(out) <- c("At low moderator", "At high moderator")
 return(out)
}


# ci.oddsratio ==============================================================
#' Confidence interval for an odds ratio
#'
#'
#' @description
#' Computes a confidence interval for an odds ratio with .5 added to each
#' cell frequency. This function requires the frequency counts from a
#' 2 x 2 contingency table for two dichotomous variables.
#'
#'
#' @param   alpha  alpha level for 1-alpha confidence
#' @param   f00    number of participants with y = 0 and x = 0
#' @param   f01    number of participants with y = 0 and x = 1
#' @param   f10    number of participants with y = 1 and x = 0
#' @param   f11    number of participants with y = 1 and x = 1
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimate of odds ratio
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @references
#' \insertRef{Fleiss2003}{statpsych}
#'
#'
#' @examples
#' ci.oddsratio(.05, 229, 28, 96, 24)
#'
#' # Should return:
#' #      Estimate        SE       LL       UL
#' # [1,] 2.044451 0.6154578 1.133267 3.688254
#'
#'
#' @importFrom stats qnorm
#' @export
ci.oddsratio <- function(alpha, f00, f01, f10, f11) {
 z <- qnorm(1 - alpha/2)
 or <- (f11 + .5)*(f00 + .5)/((f01 + .5)*(f10 + .5))
 se.lor <- sqrt(1/(f00 + .5) + 1/(f01 + .5) + 1/(f10 + .5) + 1/(f11 + .5))
 se.or <- or*se.lor
 ll <- exp(log(or) - z*se.lor)
 ul <- exp(log(or) + z*se.lor)
 out <- t(c(or, se.or, ll, ul))
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 return(out)
}


#  ci.yule ==================================================================== 
#' Confidence intervals for generalized Yule coefficients
#'
#'                            
#' @description
#' Computes confidence intervals for four generalized Yule measures of 
#' association (Yule Q, Yule Y, Digby H, and Bonett-Price Y*) using a 
#' transformation of a confidence interval for an odds ratio with .5 added to 
#' each cell frequency. This function requires the frequency counts from a 
#' 2 x 2 contingency table for two dichotomous variables. Digby H is sometimes
#' used as a crude approximation to the tetrachoric correlation. Yule Y is 
#' equal to the phi coefficient only when all marginal frequencies are equal.
#' Bonett-Price Y* is a better approximation to the phi coeffiient when the
#' marginal frequencies are not equal.
#'
#'
#' @param   alpha  alpha level for 1-alpha confidence
#' @param   f00    number of participants with y = 0 and x = 0
#' @param   f01    number of participants with y = 0 and x = 1
#' @param   f10    number of participants with y = 1 and x = 0
#' @param   f11    number of participants with y = 1 and x = 1
#'
#'
#' @references
#' \insertRef{Bonett2007}{statpsych}                   
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimate of generalized Yule coefficient
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'                                                       
#' @examples
#' ci.yule(.05, 229, 28, 96, 24)
#'
#' # Should return:
#' #      Estimate         SE         LL        UL
#' # Q:  0.3430670 0.13280379 0.06247099 0.5734020
#' # Y:  0.1769015 0.07290438 0.03126603 0.3151817
#' # H:  0.2619244 0.10514465 0.04687994 0.4537659
#' # Y*: 0.1311480 0.05457236 0.02307188 0.2361941
#'
#'
#' @importFrom stats qnorm
#' @export
ci.yule <- function(alpha, f00, f01, f10, f11) {
 z <- qnorm(1 - alpha/2)
 n <- f00 + f01 + f10 + f11
 f1 <- f00 + f01
 f2 <- f10 + f11
 f3 <- f00 + f10
 f4 <- f01 + f11
 min <- min(f1, f2, f3, f4)/n
 or <- (f11 + .5)*(f00 + .5)/((f01 + .5)*(f10 + .5))
 se.lor <- sqrt(1/(f00 + .5) + 1/(f01 + .5) + 1/(f10 + .5) + 1/(f11 + .5))
 ll.or <- exp(log(or) - z*se.lor)
 ul.or <- exp(log(or) + z*se.lor)
 Q <- (or - 1)/(or + 1)
 se.Q <- .5*(1 - Q^2)*se.lor
 ll.Q <- (ll.or - 1)/(ll.or + 1)
 ul.Q <- (ul.or - 1)/(ul.or + 1)
 Y <- (or^.5 - 1)/(or^.5 + 1)
 se.Y <- .25*(1 - Y^2)*se.lor
 ll.Y <- (ll.or^.5 - 1)/(ll.or^.5 + 1)
 ul.Y <- (ul.or^.5 - 1)/(ul.or^.5 + 1)
 H <- (or^.75 - 1)/(or^.75 + 1)
 se.H <- .375*(1 - H^2)*se.lor
 ll.H <- (ll.or^.75 - 1)/(ll.or^.75 + 1)
 ul.H <- (ul.or^.75 - 1)/(ul.or^.75 + 1)
 c <- .5 - (.5 - min)^2
 Y2 <- (or^c - 1)/(or^c + 1)
 se.Y2 <- (c/2)*(1 - Y2^2)*se.lor
 ll.Y2 <- (ll.or^c - 1)/(ll.or^c + 1)
 ul.Y2 <- (ul.or^c - 1)/(ul.or^c + 1)
 out1 <- t(c(Q, se.Q, ll.Q, ul.Q))
 out2 <- t(c(Y, se.Y, ll.Y, ul.Y))
 out3 <- t(c(H, se.H, ll.H, ul.H))
 out4 <- t(c(Y2, se.Y2, ll.Y2, ul.Y2))
 out <- rbind(out1, out2, out3, out4)
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 rownames(out) <- c("Q:", "Y:", "H:", "Y*:")
 return(out)
}
          
       
#  ci.phi ====================================================================
#' Confidence interval for a phi correlation
#'
#'
#' @description
#' Computes a confidence interval for a phi correlation. This function requires 
#' the frequency counts from a 2 x 2 contingency table for two dichotomous 
#' variables. This measure of association is usually most appropriate when both
#' dichotomous variables are naturally dichotomous.
#'
#'
#' @param   alpha  alpha level for 1-alpha confidence
#' @param   f00    number of participants with y = 0 and x = 0
#' @param   f01    number of participants with y = 0 and x = 1
#' @param   f10    number of participants with y = 1 and x = 0
#' @param   f11    number of participants with y = 1 and x = 1
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimate of phi correlation
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @references
#' \insertRef{Bishop1975}{statpsych}
#'
#'
#' @examples
#' ci.phi(.05, 229, 28, 96, 24)
#'
#' # Should return:
#' #  Estimate         SE         LL        UL
#' # 0.1229976 0.05746271 0.01037273 0.2356224
#'
#'
#' @importFrom stats qnorm
#' @export
ci.phi <- function(alpha, f00, f01, f10, f11) {
 z <- qnorm(1 - alpha/2)
 n <- f00 + f01 + f10 + f11
 p00 <- f00/n; p01 <- f01/n; p10 <- f10/n; p11 <- f11/n;
 p0x <- (f00 + f01)/n; p1x <- (f10 + f11)/n
 px0 <- (f00 + f10)/n; px1 <- (f01 + f11)/n
 phi <- (p11*p00 - p10*p01)/sqrt(p1x*p0x*px1*px0)
 v1 <- 1 - phi^2 
 v2 <- phi + .5*phi^3
 v3 <- (p0x - p1x)*(px0 - px1)/sqrt(p0x*p1x*px0*px1) 
 v4 <- (.75*phi^2)*((p0x - p1x)^2/(p0x*p1x) + (px0 - px1)^2/(px0*px1))
 se <- sqrt((v1 + v2*v3 + v4)/n)
 ll <- phi - z*se
 ul <- phi + z*se
 out <- t(c(phi, se, ll, ul))
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 rownames(out) <- ""
 return(out)
}


#  ci.biphi ==================================================================
#' Confidence interval for a biserial-phi correlation
#'
#'
#' @description
#' Computes a confidence interval for a biserial-phi correlation using a
#' transformation of a confidence interval for an odds ratio with .5 added to
#' each cell frequency. This measure of association assumes the group variable
#' is naturally dichotomous and the response variable is artificially
#' dichotomous. 
#'
#'
#' @param   alpha  alpha level for 1-alpha confidence
#' @param   f1     number of participants in group 1 who have the attribute
#' @param   f2     number of participants in group 2 who have the attribute
#' @param   n1     sample size for group 1
#' @param   n2     sample size for group 2
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimate of biserial-phi correlation
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @references
#' \insertRef{Ulrich2004}{statpsych}
#'
#'
#' @examples
#' ci.biphi(.05, 46, 15, 100, 100)
#'
#' # Should return:
#' #  Estimate         SE        LL       UL
#' # 0.4145733 0.07551281 0.2508866 0.546141
#'
#'
#' @importFrom stats qnorm
#' @export
ci.biphi <- function(alpha, f1, f2, n1, n2) {
 if (f1 > n1) {stop("f cannot be greater than n")}
 if (f2 > n2) {stop("f cannot be greater than n")}
 z <- qnorm(1 - alpha/2)
 f00 <- f1
 f10 <- n1 - f1
 f01 <- f2
 f11 <- n2 - f2
 p1 <- n1/(n1 + n2)
 p2 <- n2/(n1 + n2)
 or <- (f11 + .5)*(f00 + .5)/((f01 + .5)*(f10 + .5))
 lor <- log(or)
 se.lor <- sqrt(1/(f00 + .5) + 1/(f01 + .5) + 1/(f10 + .5) + 1/(f11 + .5))
 LL1 <- lor - z*se.lor
 UL1 <- lor + z*se.lor
 c <- 2.89/(p1*p2)
 biphi <- lor/sqrt(lor^2 + c)
 se.biphi <- sqrt(c^2/(lor^2 + c)^3)*se.lor
 ll <- LL1/sqrt(LL1^2 + 2.89/(p1*p2))
 ul <- UL1/sqrt(UL1^2 + 2.89/(p1*p2))
 out <- t(c(biphi, se.biphi, ll, ul))
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 rownames(out) <- ""
 return(out)
}


#  ci.tetra ==================================================================
#' Confidence interval for a tetrachoric correlation 
#'
#'
#' @description
#' Computes a confidence interval for an approximation to the tetrachoric
#' correlation. This function requires the frequency counts from a 2 x 2 
#' contingency table for two dichotomous variables. This measure of 
#' association assumes both of the dichotomous variables are artificially 
#' dichotomous. An approximate standard error is recovered from the
#' confidence interval.
#'
#'
#' @param   alpha  alpha level for 1-alpha confidence
#' @param   f00    number of participants with y = 0 and x = 0
#' @param   f01    number of participants with y = 0 and x = 1
#' @param   f10    number of participants with y = 1 and x = 0
#' @param   f11    number of participants with y = 1 and x = 1
#'
#'
#' @references
#' \insertRef{Bonett2005}{statpsych}
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimate of tetrachoric approximation
#' * SE - recovered standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @examples
#' ci.tetra(.05, 46, 15, 54, 85)
#'
#' # Should return:
#' #  Estimate         SE        LL        UL
#' # 0.5135167 0.09301703 0.3102345 0.6748546
#'
#'
#' @importFrom stats qnorm
#' @export
ci.tetra <- function(alpha, f00, f01, f10, f11) {
 z <- qnorm(1 - alpha/2)
 n <- f00 + f01 + f10 + f11
 or <- (f11 + .5)*(f00 + .5)/((f01 + .5)*(f10 + .5))
 r1 <- (f00 + f01 + 1)/(n + 2)
 r2 <- (f10 + f11 + 1)/(n + 2)
 c1 <- (f00 + f10 + 1)/(n + 2)
 c2 <- (f01 + f11 + 1)/(n + 2)
 pmin <- min(c1, c2, r1, r2)
 c <- (1 - abs(r1 - c1)/5 - (.5 - pmin)^2)/2
 lor <- log(or)
 se.lor <- sqrt(1/(f00 + .5) + 1/(f01 + .5) + 1/(f10 + .5) + 1/(f11 + .5))
 LL1 <- exp(lor - z*se.lor)
 UL1 <- exp(lor + z*se.lor)
 tetra <- cos(3.14159/(1 + or^c))
 ll <- cos(3.14159/(1 + LL1^c))
 ul <- cos(3.14159/(1 + UL1^c))
 se <- (ul - ll)/(2*z)
 out <- t(c(tetra, se, ll, ul))
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 rownames(out) <- ""
 return(out)
}


#  ci.kappa ================================================================== 
#' Confidence interval for two kappa reliability coefficients
#'
#'
#' @description
#' Computes confidence intervals for the intraclass kappa coefficient and
#' Cohen's kappa coefficient with two dichotomous ratings. 
#'
#'
#' @param   alpha  alpha level for 1-alpha confidence
#' @param   f00    number of objects rated y = 0 and x = 0
#' @param   f01    number of objects rated y = 0 and x = 1
#' @param   f10    number of objects rated y = 1 and x = 0
#' @param   f11    number of objects rated y = 1 and x = 1
#'
#'
#' @return
#' Returns a 2-row matrix. The results in row 1 are for the intraclass
#' kappa. The results in row 2 are for Cohen's kappa. The columns are:
#' * Estimate - estimate of interrater reliability 
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @references
#' \insertRef{Fleiss2003}{statpsych}
#'
#'
#' @examples
#' ci.kappa(.05, 31, 12, 4, 58)
#'
#' # Should return:
#' #               Estimate         SE        LL        UL
#' # IC kappa:    0.6736597 0.07479965 0.5270551 0.8202643
#' # Cohen kappa: 0.6756757 0.07344761 0.5317210 0.8196303
#'
#'
#' @importFrom stats qnorm
#' @export
ci.kappa <- function(alpha, f00, f01, f10, f11) {
 z <- qnorm(1 - alpha/2)
 n <- f00 + f01 + f10 + f11 
 p00 <- f00/n; p01 <- f01/n; p10 <- f10/n; p11 <- f11/n
 p1 <- (2*f00 + f01 + f10)/(2*n)
 k1 <- 4*(p00*p11 - p01*p10) - (p01 - p10)^2
 k2 <- (2*p00 + p01 + p10)*(2*p11 + p01 + p10)
 k <- k1/k2 
 se.k <- sqrt(((1 - k)/n)*((1 - k)*(1 - 2*k) + k*(2 - k)/(2*p1*(1 - p1))))
 LL.k <- k - z*se.k
 UL.k <- k + z*se.k 
 pr <- (p00 + p01)*(p00 + p10) + (p10 + p11)*(p01 + p11)
 c <- ((p00 + p11) - pr)/(1 - pr)
 a1 <- p11*(1 - (p10 + p01 + 2*p11)*(1 - c))^2 + p00*(1 - (p10 + p01 + 2*p00)*(1 - c))^2
 a2 <- p10*(p11 + p00 + 2*p01)^2*(1 - c)^2 + p01*(p11 + p00 + 2*p10)^2*(1 - c)^2
 a3 <- (c - pr*(1 - c))^2
 se.c <- sqrt(a1 + a2 - a3)/((1 - pr)*sqrt(n))
 LL.c <- c - z*se.c
 UL.c <- c + z*se.c 
 out1 <- c(k, se.k, LL.k, UL.k)
 out2 <- c(c, se.c, LL.c, UL.c)
 out <- rbind(out1, out2)
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 rownames(out) <- c("IC kappa:", "Cohen kappa:")
 return(out)
}


#  ci.agree ==================================================================
#' Confidence interval for a G-index of agreement
#'
#'
#' @description
#' Computes an adjusted Wald  confidence interval for a G-index of agreement
#' between two polychotomous ratings. This function requires the number of 
#' objects that were given the same rating by both raters. The G-index
#' corrects for chance agreement.
#'
#'
#' @param   alpha  alpha level for 1-alpha confidence
#' @param   n      sample size
#' @param   f      number of objects rated in agreement
#' @param   k      number of rating categories
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - maximum likelihood estimate of G-index 
#' * SE - standard error
#' * LL - lower limit of the adjusted Wald confidence interval
#' * UL - upper limit of the adjusted Wald confidence interval
#'
#'
#' @references
#' \insertRef{Bonett2022}{statpsych}
#'
#'
#' @examples
#' ci.agree(.05, 100, 80, 4)
#'
#' # Should return:
#' #  Estimate         SE        LL        UL
#' # 0.7333333 0.05333333 0.6132949 0.8226025
#'
#'
#' @importFrom stats qnorm
#' @export
ci.agree <- function(alpha, n, f, k) {
 if (f > n) {stop("f cannot be greater than n")}
 z <- qnorm(1 - alpha/2)
 a <- k/(k - 1)
 g.mle <- a*f/n - 1/(k - 1)
 p.mle <- f/n
 p.adj <- (f + 2)/(n + 4)
 se.g <- a*sqrt(p.mle*(1 - p.mle)/n)
 LL.g <- a*(p.adj - z*sqrt(p.adj*(1 - p.adj)/(n + 4))) - 1/(k - 1) 
 UL.g <- a*(p.adj + z*sqrt(p.adj*(1 - p.adj)/(n + 4))) - 1/(k - 1) 
 out <- t(c(g.mle, se.g, LL.g, UL.g))
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 rownames(out) <- ""
 return(out)
}


#  ci.agree2 =================================================================
#' Confidence interval for G-index difference in a 2-group design
#'
#'                          
#' @description
#' Computes adjusted Wald confidence intervals for the G-index of agreement 
#' within each group and the difference of G-indices. 
#'
#'
#' @param  alpha   alpha level for simultaneous 1-alpha confidence
#' @param  n1      sample size (objects) in group 1
#' @param  f1      number of objects rated in agreement in group 1
#' @param  n2      sample size (objects) in group 2
#' @param  f2      number of objects rated in agreement in group 2
#' @param  r       number of rating categories
#'
#'
#' @return
#' Returns a 3-row matrix. The rows are:
#' * Row 1: G-index for group 1
#' * Row 2: G-index for group 2
#' * Row 3: G-index difference
#'
#'
#' The columns are:
#' * Estimate - maximum likelihood estimate of G-index and difference  
#' * SE - standard error
#' * LL - lower limit of adjusted Wald confidence interval
#' * UL - upper limit of adjusted Wald confidence interval
#'
#'
#' @references
#' \insertRef{Bonett2022}{statpsych}
#'
#'
#' @examples
#' ci.agree2(.05, 75, 70, 60, 45, 2)
#'
#' # Should return:
#' #          Estimate         SE        LL        UL
#' # G1      0.8666667 0.02880329 0.6974555 0.9481141
#' # G2      0.5000000 0.05590170 0.2523379 0.6851621
#' # G1 - G2 0.3666667 0.06288585 0.1117076 0.6088621
#'
#'
#' @importFrom stats qnorm
#' @export
ci.agree2 <- function(alpha, n1, f1, n2, f2, r) {
 if (f1 > n1) {stop("f cannot be greater than n")}
 if (f2 > n2) {stop("f cannot be greater than n")}
 z <- qnorm(1 - alpha/2)
 a <- r/(r - 1)
 p1.ml <- f1/n1
 p1 <- (f1 + 2)/(n1 + 4)
 G1 <- a*p1.ml - 1/(r - 1)
 se1 <- sqrt(p1*(1 - p1)/(n1 + 4))
 se1.ml <- sqrt(p1.ml*(1 - p1.ml)/n1)
 LL1 <- a*(p1 - z*se1) - 1/(r - 1)
 UL1 <- a*(p1 + z*se1) - 1/(r - 1) 
 p2.ml <- f2/n2
 p2 <- (f2 + 2)/(n2 + 4)
 G2 <- a*p2.ml - 1/(r - 1)
 se2 <- sqrt(p2*(1 - p2)/(n2 + 4))
 se2.ml <- sqrt(p2.ml*(1 - p2.ml)/n2)
 LL2 <- a*(p2 - z*se2) - 1/(r - 1)
 UL2 <- a*(p2 + z*se2) - 1/(r - 1) 
 p1.d <- (f1 + 1)/(n1 + 2)
 p2.d <- (f2 + 1)/(n2 + 2)
 se.d <- sqrt(p1.d*(1 - p1.d)/(n1 + 2) + p2.d*(1 - p2.d)/(n2 + 2))
 se.d.ml <- sqrt(se1.ml^2 + se2.ml^2)
 LL3 <- a*(p1.d - p2.d - z*se.d)
 UL3 <- a*(p1.d - p2.d + z*se.d) 
 out1 <- t(c(G1, se1.ml, LL1, UL1))
 out2 <- t(c(G2, se2.ml, LL2, UL2))
 out3 <- t(c(G1 - G2, se.d.ml, LL3, UL3))
 out <- rbind(out1, out2, out3)
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 rownames(out) <- c("G1", "G2", "G1 - G2")
 return(out)
}
          

# ci.agree.3rater =============================================================
#' Computes confidence intervals for a 3-rater design with dichotomous ratings
#'
#'     
#' @description
#' Computes adjusted Wald confidence intervals for a G-index of agreement for 
#' all pairs of raters in a 3-rater design with a dichotomous rating, and 
#' computes adjusted Wald confidence intervals for differences of all pairs of 
#' G agreement. An adjusted Wald confidence interval for unanimous G agreement 
#' among the three raters is also computed. In the three-rater design, 
#' unanimous G agreement is equal to the average of all pairs of G agreement. 
#' The G-index corrects for chance agreement.
#'
#'  
#' @param  alpha    alpha level for 1-alpha confidence
#' @param  f        vector of frequency counts from 2x2x2 table where
#'                  f = \[ f111, f112, f121, f122, f211, f212, f221, f222 \],
#'                  first subscript represents the rating of rater 1,
#'                  second subscript represents the rating of rater 2, and
#'                  third subscript represents the rating of rater 3
#'
#' 
#' @references
#' \insertRef{Bonett2022}{statpsych}
#'
#'
#' @return 
#' Returns a 7-row matrix. The rows are:
#' * G(1,2): G-index for raters 1 and 2
#' * G(1,3): G-index for raters 1 and 3
#' * G(2,3): G-index for raters 2 and 3
#' * G(1,2)-G(1,3): difference in G(1,2) and G(1,3)
#' * G(1,2)-G(2,3): difference in G(1,2) and G(2,3)
#' * G(2,3)-G(1,3): difference in G(2,3) and G(1,3)
#' * G(3): G-index of unanimous agreement for all three raters
#'
#'
#' The columns are:
#' * Estimate - estimate of G-index (two-rater, difference, or unanimous)  
#' * LL - lower limit of adjusted Wald confidence interval
#' * UL - upper limit of adjusted Wald confidence interval
#'
#' 
#' 
#' @examples
#' f <- c(100, 6, 4, 40, 20, 1, 9, 120)
#' ci.agree.3rater(.05, f)
#'
#' # Should return:
#' #                  Estimate          LL         UL
#' # G(1,2)         0.56666667  0.46601839  0.6524027
#' # G(1,3)         0.50000000  0.39564646  0.5911956
#' # G(2,3)         0.86666667  0.79701213  0.9135142
#' # G(1,2)-G(1,3)  0.06666667  0.00580397  0.1266464
#' # G(1,2)-G(2,3) -0.30000000 -0.40683919 -0.1891873
#' # G(2,3)-G(1,3) -0.36666667 -0.46222023 -0.2662566
#' # G(3)           0.64444444  0.57382971  0.7068720
#'  
#' 
#' @importFrom stats qnorm
#' @export         
ci.agree.3rater <- function(alpha, f) {
 z <- qnorm(1 - alpha/2)
 n <- sum(f)
 f111 <- f[1]
 f112 <- f[2]
 f121 <- f[3]
 f122 <- f[4]
 f211 <- f[5]
 f212 <- f[6]
 f221 <- f[7]
 f222 <- f[8]
 p12.ml <- (f111 + f112 + f221 + f222)/n 
 p12 <- (f111 + f112 + f221 + f222 + 2)/(n + 4)
 p13.ml <- (f111 + f121 + f212 + f222)/n
 p13 <- (f111 + f121 + f212 + f222 + 2)/(n + 4)
 p23.ml <- (f111 + f211 + f122 + f222)/n
 p23 <- (f111 + f211 + f122 + f222 + 2)/(n + 4)
 G12.ml <- 2*p12.ml - 1
 G12 <- 2*p12 - 1
 G13.ml <- 2*p13.ml - 1
 G13 <- 2*p13 - 1
 G23.ml <- 2*p23.ml - 1
 G23 <- 2*p23 - 1
 se.G12 <- sqrt(p12*(1 - p12)/(n + 4))
 se.G13 <- sqrt(p13*(1 - p13)/(n + 4))
 se.G23 <- sqrt(p23*(1 - p23)/(n + 4))
 p1.ml <- (f112 + f221)/n
 p1 <- (f112 + f221 + 1)/(n + 2)
 p2.ml <- (f121 + f212)/n
 p2 <- (f121 + f212 + 1)/(n + 2)
 p3.ml <- (f211 + f122)/n
 p3 <- (f211 + f122 + 1)/(n + 2)
 G12_13.ml <- 2*(p1.ml - p2.ml)
 G12_13 <- 2*(p1 - p2)
 G12_23.ml <- 2*(p1.ml - p3.ml)
 G12_23 <- 2*(p1 - p3)
 G13_23.ml <- 2*(p2.ml - p3.ml)
 G13_23 <- 2*(p2 - p3)
 se.G12_13 <- sqrt((p1 + p2 - (p1 - p2)^2)/(n + 2))
 se.G12_23 <- sqrt((p1 + p3 - (p1 - p3)^2)/(n + 2))
 se.G13_23 <- sqrt((p2 + p3 - (p2 - p3)^2)/(n + 2))
 p123.ml <- (f111 + f222)/n
 p123 <- (f111 + f222 + 2)/(n + 4)
 G3.ml <- (4*p123.ml - 1)/3
 G3 <- (4*p123 - 1)/3
 se.G3 <- sqrt(p123*(1 - p123)/(n + 4))
 LL.G12 <- 2*(p12 - z*se.G12) - 1
 UL.G12 <- 2*(p12 + z*se.G12) - 1
 LL.G13 <- 2*(p13 - z*se.G13) - 1
 UL.G13 <- 2*(p13 + z*se.G13) - 1
 LL.G23 <- 2*(p23 - z*se.G23) - 1
 UL.G23 <- 2*(p23 + z*se.G23) - 1
 LL.G12_13 <- 2*(p1 - p2 - z*se.G12_13)
 UL.G12_13 <- 2*(p1 - p2 + z*se.G12_13)
 LL.G12_23 <- 2*(p1 - p3 - z*se.G12_23)
 UL.G12_23 <- 2*(p1 - p3 + z*se.G12_23)
 LL.G13_23 <- 2*(p2 - p3 - z*se.G13_23)
 UL.G13_23 <- 2*(p2 - p3 + z*se.G13_23)
 LL.G3 <- (4/3)*(p123 - z*se.G3) - 1/3
 UL.G3 <- (4/3)*(p123 + z*se.G3) - 1/3
 out1 <- t(c(G12.ml, LL.G12, UL.G12))
 out2 <- t(c(G13.ml, LL.G13, UL.G13))
 out3 <- t(c(G23.ml, LL.G23, UL.G23))
 out4 <- t(c(G12_13.ml, LL.G12_13, UL.G12_13))
 out5 <- t(c(G12_23.ml, LL.G12_23, UL.G12_23))
 out6 <- t(c(G13_23.ml, LL.G13_23, UL.G13_23))
 out7 <- t(c(G3.ml, LL.G3, UL.G3))
 out <- rbind(out1, out2, out3, out4, out5, out6, out7)
 colnames(out) <- c("Estimate", "LL", "UL")
 rownames(out) <-  c("G(1,2)", "G(1,3)", "G(2,3)", "G(1,2)-G(1,3)", "G(1,2)-G(2,3)",
                     "G(2,3)-G(1,3)","G(3)")
 return(out)
}

          
# ci.popsize ================================================================= 
#' Confidence interval for an unknown population size
#'
#'
#' @description
#' Computes a Wald confidence interval for an unknown population size using 
#' mark-recapture sampling. This method assumes independence of the two
#' samples. This function requires the frequency counts from an incomplete
#' 2 x 2 contingency table for the two samples (f11 is the unknown number
#' of people who were not observed in either sample). This method sets the
#' estimated odds ratio (with .5 added to each cell) to 1 and solves for
#' unobserved cell frequency. An approximate standard error is recovered
#' from the confidence interval.
#'
#'
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  f00    number of people observed in both samples
#' @param  f01    number of people observed in first sample but not second sample
#' @param  f10    number of people observed in second sample but not first sample
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimate of the unknown population size 
#' * SE - recovered standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @examples
#' ci.popsize(.05, 794, 710, 741)
#'
#' # Should return:
#' # Estimate       SE   LL   UL
#' #     2908 49.49071 2818 3012
#'
#'
#' @importFrom stats qnorm
#' @export
ci.popsize <- function(alpha, f00, f01, f10) {
 z <- qnorm(1 - alpha/2)
 n0 <- f00 + f01 + f10
 f11 <- (f01 + .5)*(f10 + .5)/(f00 + .5) - .5
 se <- sqrt(1/(f00 + .5) + 1/(f01 + .5) + 1/(f10 + .5) + 1/f11)
 N <- round(n0 + f11)
 ll <- round(n0 + exp(log(f11) - z*se))
 ul <- round(n0 + exp(log(f11) + z*se))
 se <- (ul - ll)/(2*z)
 out <- t(c(N, se, ll, ul))
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 rownames(out) <- ""
 return(out)
}


#  ci.cramer  ======================================================================
#' Confidence interval for Cramer's V
#'
#'
#' @description
#' Computes a confidence interval for a population Cramer's V coefficient
#' of nominal association for an r x s contingency table and its approximate
#' standard error. The confidence interval is based on a noncentral chi-square 
#' distribution, and an approximate standard error is recovered from the
#' confidence interval.
#'
#'
#' @param  alpha    alpha value for 1-alpha confidence
#' @param  chisqr   Pearson chi-square test statistic of independence
#' @param  r        number of rows in contingency table
#' @param  c        number of columns in contingency table
#' @param  n        sample size
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimate of Cramer's V 
#' * SE - recovered standard error 
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @references
#' \insertRef{Smithson2003}{statpsych}
#'
#'
#' @examples
#' ci.cramer(.05, 19.21, 2, 3, 200)
#'
#' # Should return:
#' # Estimate     SE     LL     UL
#' #   0.3099 0.0718 0.1601 0.4417
#'  
#' 
#' @importFrom stats pchisq
#' @importFrom stats qnorm
#' @export
ci.cramer <- function(alpha, chisqr, r, c, n) {
 alpha1 <- alpha/2
 alpha2 <- 1 - alpha/2
 z <- qnorm(1 - alpha/2)
 k <- min(r - 1, c - 1)
 df <- (r - 1)*(c - 1)
 v <- sqrt(chisqr/(n*k))
 du <- n*k 
 nc <- seq(0, du, by = .001)
 p <- pchisq(chisqr, df, nc)
 k1 <- which(min(abs(p - alpha2)) == abs(p - alpha2))[[1]]
 dL <- nc[k1]
 # version 1.5 corrects an error in the Smithson CI equation 4.13 
 ll <- sqrt(dL/(n*k))
 k2 <- which(min(abs(p - alpha1)) == abs(p - alpha1))[[1]]
 dU <- nc[k2]
 ul <- sqrt(dU/(n*k))
 se <- (ul - ll)/(2*z)
 out <- round(t(c(v, se, ll, ul)), 4)
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 rownames(out) <- ""
 return(out)
}


#  ci.2x2.prop.bs ==========================================================
#' Computes tests and confidence intervals of effects in a 2x2 between-
#' subjects design for proportions 
#'
#'
#' @description
#' Computes adjusted Wald confidence intervals and tests for the AB 
#' interaction effect, main effect of A, main efect of B, simple main effects
#' of A, and simple main effects of B in a 2x2 between-subjects factorial 
#' design with a dichotomous response variable. The input vector of 
#' frequency counts is f = \[ f11, f12, f21, f22 \], and the input vector of 
#' sample sizes is n = \[ n11, n12, n21, n22 \] where the first subscript 
#' represents the levels of Factor A and the second subscript represents the 
#' levels of Factor B.
#'
#'
#' @param   alpha   alpha level for 1-alpha confidence
#' @param   f       vector of frequency counts of participants who have the attribute
#' @param   n       vector of sample sizes
#'
#'
#' @return
#' Returns a 7-row matrix (one row per effect). The columns are:
#' * Estimate - adjusted estimate of effect
#' * SE - standard error 
#' * z - z test statistic for test of null hypothesis
#' * p - two-sided p-value 
#' * LL - lower limit of the adjusted Wald confidence interval
#' * UL - upper limit of the adjusted Wald confidence interval
#'
#'
#' @references
#' \insertRef{Price2004}{statpsych}
#'
#'
#' @examples
#' f <- c(15, 24, 28, 23)
#' n <- c(50, 50, 50, 50)
#' ci.2x2.prop.bs(.05, f, n)
#'
#' # Should return:
#' #             Estimate         SE          z           p          LL          UL
#' # AB:      -0.27450980 0.13692496 -2.0048193 0.044982370 -0.54287780 -0.00614181
#' # A:       -0.11764706 0.06846248 -1.7184165 0.085720668 -0.25183106  0.01653694
#' # B:       -0.03921569 0.06846248 -0.5728055 0.566776388 -0.17339968  0.09496831
#' # A at b1: -0.25000000 0.09402223 -2.6589456 0.007838561 -0.43428019 -0.06571981
#' # A at b2:  0.01923077 0.09787658  0.1964798 0.844234654 -0.17260380  0.21106534
#' # B at a1: -0.17307692 0.09432431 -1.8349132 0.066518551 -0.35794917  0.01179533
#' # B at a2:  0.09615385 0.09758550  0.9853293 0.324462356 -0.09511021  0.28741790
#'
#'
#' @importFrom stats qnorm
#' @importFrom stats pnorm
#' @export
ci.2x2.prop.bs <- function(alpha, f, n) {
 s <- sum(as.integer(as.logical(n < f)))
 if (s > 0) {stop("f cannot be greater than n")}
 zcrit <- qnorm(1 - alpha/2)
 v1 <- c(1, -1, -1, 1)
 v2 <- c(.5, .5, -.5, -.5)
 v3 <- c(.5, -.5, .5, -.5)
 v4 <- c(1, 0, -1, 0)
 v5 <- c(0, 1, 0, -1)
 v6 <- c(1, -1, 0, 0)
 v7 <- c(0, 0, 1, -1)
 p.4 <- (f + .5)/(n + 1)
 p.2 <- (f + 1)/(n + 2)
 est1 <- t(v1)%*%p.4
 se1 <- sqrt(t(v1)%*%diag(p.4*(1 - p.4))%*%solve(diag(n + 1))%*%v1)
 z1 <- est1/se1
 p1 <- 2*(1 - pnorm(abs(z1)))
 LL1 <- est1 - zcrit*se1
 UL1 <- est1 + zcrit*se1
 row1 <- c(est1, se1, z1, p1, LL1, UL1)
 est2 <- t(v2)%*%p.4
 se2 <- sqrt(t(v2)%*%diag(p.4*(1 - p.4))%*%solve(diag(n + 1))%*%v2)
 z2 <- est2/se2
 p2 <- 2*(1 - pnorm(abs(z2)))
 LL2 <- est2 - zcrit*se2
 UL2 <- est2 + zcrit*se2
 row2 <- c(est2, se2, z2, p2, LL2, UL2)
 est3 <- t(v3)%*%p.4
 se3 <- sqrt(t(v3)%*%diag(p.4*(1 - p.4))%*%solve(diag(n + 1))%*%v3)
 z3 <- est3/se3
 p3 <- 2*(1 - pnorm(abs(z3)))
 LL3 <- est3 - zcrit*se3
 UL3 <- est3 + zcrit*se3
 row3 <- c(est3, se3, z3, p3, LL3, UL3)
 est4 <- t(v4)%*%p.2
 se4 <- sqrt(t(v4)%*%diag(p.2*(1 - p.2))%*%solve(diag(n + 2))%*%v4)
 z4 <- est4/se4
 p4 <- 2*(1 - pnorm(abs(z4)))
 LL4 <- est4 - zcrit*se4
 UL4 <- est4 + zcrit*se4
 row4 <- c(est4, se4, z4, p4, LL4, UL4)
 est5 <- t(v5)%*%p.2
 se5 <- sqrt(t(v5)%*%diag(p.2*(1 - p.2))%*%solve(diag(n + 2))%*%v5)
 z5 <- est5/se5
 p5 <- 2*(1 - pnorm(abs(z5)))
 LL5 <- est5 - zcrit*se5
 UL5 <- est5 + zcrit*se5
 row5 <- c(est5, se5, z5, p5, LL5, UL5)
 est6 <- t(v6)%*%p.2
 se6 <- sqrt(t(v6)%*%diag(p.2*(1 - p.2))%*%solve(diag(n + 2))%*%v6)
 z6 <- est6/se6
 p6 <- 2*(1 - pnorm(abs(z6)))
 LL6 <- est6 - zcrit*se6
 UL6 <- est6 + zcrit*se6
 row6 <- c(est6, se6, z6, p6, LL6, UL6)
 est7 <- t(v7)%*%p.2
 se7 <- sqrt(t(v7)%*%diag(p.2*(1 - p.2))%*%solve(diag(n + 2))%*%v7)
 z7 <- est7/se7
 p7 <- 2*(1 - pnorm(abs(z7)))
 LL7 <- est7 - zcrit*se7
 UL7 <- est7 + zcrit*se7
 row7 <- c(est7, se7, z7, p7, LL7, UL7)
 out <- rbind(row1, row2, row3, row4, row5, row6, row7)
 rownames(out) <- c("AB:", "A:", "B:", "A at b1:", "A at b2:", "B at a1:", "B at a2:")
 colnames(out) = c("Estimate", "SE", "z", "p", "LL", "UL")
 return(out)
}


#  ci.2x2.prop.mixed ==========================================================
#' Computes tests and confidence intervals of effects in a 2x2 mixed factorial
#' design for proportions 
#'
#'
#' @description
#' Computes adjusted Wald confidence intervals and tests for the AB 
#' interaction effect, main effect of A, main efect of B, simple main effects
#' of A, and simple main effects of B in a 2x2 mixed factorial design with a
#' dichotomous response variable where Factor A is a within-subjects factor 
#' and Factor B is a between-subjects factor. The 4x1 vector of frequency 
#' counts for Factor A within each group is f00, f01, f10, f11 where fij is 
#' the number of participants with a response of i = 0 or 1 at level 1 of 
#' Factor A and a response of j = 0 or 1 at level 2 of Factor A. 
#'
#'
#' @param   alpha   alpha level for 1-alpha confidence
#' @param   group1  vector of frequency counts from 2x2 contingency table for Factor A in group 1
#' @param   group2  vector of frequency counts from 2x2 contingency table for Factor A in group 2
#'
#'
#' @return
#' Returns a 7-row matrix (one row per effect). The columns are:
#' * Estimate - adjusted estimate of effect
#' * SE - standard error of estimate
#' * z - z test statistic 
#' * p - two-sided p-value
#' * LL - lower limit of the adjusted Wald confidence interval
#' * UL - upper limit of the adjusted Wald confidence interval
#'
#'
#' @examples
#' group1 <- c(125, 14, 10, 254)
#' group2 <- c(100, 16, 9, 275)
#' ci.2x2.prop.mixed (.05, group1, group2)
#'
#' # Should return:
#' #              Estimate          SE          z          p          LL           UL
#' # AB:       0.007555369 0.017716073  0.4264697 0.66976559 -0.02716750  0.042278234
#' # A:       -0.013678675 0.008858036 -1.5442107 0.12253730 -0.03104011  0.003682758
#' # B:       -0.058393219 0.023032656 -2.5352360 0.01123716 -0.10353640 -0.013250043
#' # A at b1: -0.009876543 0.012580603 -0.7850612 0.43241768 -0.03453407  0.014780985
#' # A at b2: -0.017412935 0.012896543 -1.3502018 0.17695126 -0.04268969  0.007863824
#' # B at a1: -0.054634236 0.032737738 -1.6688458 0.09514794 -0.11879902  0.009530550
#' # B at a2: -0.062170628 0.032328556 -1.9230871 0.05446912 -0.12553343  0.001192177
#'
#'
#' @importFrom stats qnorm
#' @importFrom stats pnorm
#' @export
ci.2x2.prop.mixed <- function(alpha, group1, group2) {
 zcrit <- qnorm(1 - alpha/2)
 n1 <- sum(group1)
 n2 <- sum(group2)
 f11 <- group1[1]; f12 <- group1[2]; f13 <- group1[3]; f14 <- group1[4]
 f21 <- group2[1]; f22 <- group2[2]; f23 <- group2[3]; f24 <- group2[4]
 p1 <- (f12 + .5)/(n1 + 1)
 p2 <- (f13 + .5)/(n1 + 1)
 p3 <- (f22 + .5)/(n2 + 1)
 p4 <- (f23 + .5)/(n2 + 1)
 est1 <- (p2 - p1) - (p4 - p3)
 v1 <- (p1 + p2 - (p1 - p2)^2)/(n1 + 1)
 v2 <- (p3 + p4 - (p3 - p4)^2)/(n2 + 1)
 se1 <- sqrt(v1 + v2)
 z1 <- est1/se1
 pval1 <- 2*(1 - pnorm(abs(z1)))
 LL1 <- est1 - zcrit*se1
 UL1 <- est1 + zcrit*se1
 row1 <- c(est1, se1, z1, pval1, LL1, UL1)
 est2 <- ((p2 - p1) + (p4 - p3))/2
 se2 <- se1/2
 z2 <- est2/se2
 pval2 <- 2*(1 - pnorm(abs(z2)))
 LL2 <- est2 - zcrit*se2
 UL2 <- est2 + zcrit*se2
 row2 <- c(est2, se2, z2, pval2, LL2, UL2)
 p1 <- (2*f14 + f12 + f13 + 1)/(2*(n1 + 2))
 p2 <- (2*f24 + f22 + f23 + 1)/(2*(n2 + 2))
 est3 <- p1 - p2
 se3 <- sqrt(p1*(1 - p1)/(2*(n1 + 2)) + p2*(1 - p2)/(2*(n2 + 2)))
 z3 <- est3/se3
 pval3 <- 2*(1 - pnorm(abs(z3)))
 LL3 <- est3 - zcrit*se3
 UL3 <- est3 + zcrit*se3
 row3 <- c(est3, se3, z3, pval3, LL3, UL3)
 p1 <- (f12 + 1)/(n1 + 2)
 p2 <- (f13 + 1)/(n1 + 2)
 est4 <- p2 - p1
 se4 <- sqrt((p1 + p2 - (p1 - p2)^2)/(n1 + 2))
 z4 <- est4/se4
 pval4 <- 2*(1 - pnorm(abs(z4)))
 LL4 <- est4 - zcrit*se4
 UL4 <- est4 + zcrit*se4
 row4 <- c(est4, se4, z4, pval4, LL4, UL4)
 p1 <- (f22 + 1)/(n2 + 2)
 p2 <- (f23 + 1)/(n2 + 2)
 est5 <- p2 - p1
 se5 <- sqrt((p1 + p2 - (p1 - p2)^2)/(n2 + 2))
 z5 <- est5/se5
 pval5 <- 2*(1 - pnorm(abs(z5)))
 LL5 <- est5 - zcrit*se5
 UL5 <- est5 + zcrit*se5
 row5 <- c(est5, se5, z5, pval5, LL5, UL5)
 p1 <- (f14 + f13 + 1)/(n1 + 2)
 p2 <- (f24 + f23 + 1)/(n2 + 2)
 est6 <- p1 - p2
 se6 <- sqrt(p1*(1 - p1)/(n1 + 2) + p2*(1 - p2)/(n2 + 2))
 z6 <- est6/se6
 pval6 <- 2*(1 - pnorm(abs(z6)))
 LL6 <- est6 - zcrit*se6
 UL6 <- est6 + zcrit*se6
 row6 <- c(est6, se6, z6, pval6, LL6, UL6)
 p1 <- (f14 + f12 + 1)/(n1 + 2)
 p2 <- (f24 + f22 + 1)/(n2 + 2)
 est7 <- p1 - p2
 se7 <- sqrt(p1*(1 - p1)/(n1 + 2) + p2*(1 - p2)/(n2 + 2))
 z7 <- est7/se7
 pval7 <- 2*(1 - pnorm(abs(z7)))
 LL7 <- est7 - zcrit*se7
 UL7 <- est7 + zcrit*se7
 row7 <- c(est7, se7, z7, pval7, LL7, UL7)
 out <- rbind(row1, row2, row3, row4, row5, row6, row7)
 rownames(out) <- c("AB:", "A:", "B:", "A at b1:", "A at b2:", "B at a1:", "B at a2:")
 colnames(out) = c("Estimate", "SE", "z", "p", "LL", "UL")
 return(out)
}


# ci.bayes.prop ==============================================================
#' Bayesian credible interval for a single proportion
#'
#'
#' @description
#' Computes a Bayesian credible interval for a single proportion using the 
#' mean and standard deviation of a prior Beta distribution along with sample
#' information. The mean and standard deviation of the posterior Beta 
#' distribution are also reported. For a noninformative prior, set the prior 
#' mean to .5 and the prior standard deviation to 1/sqrt(12) (which 
#' corresponds to a Beta(1,1) distribution). The prior variance must be 
#' less than m(1 - m) where m is the prior mean.
#'
#'
#' @param   alpha        alpha level for 1-alpha credibility interval
#' @param   prior.mean   mean of prior Beta distribution    
#' @param   prior.sd     standard deviation of prior Beta distribution 
#' @param   f            number of participants who have the attribute
#' @param   n            sample size
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Posterior mean - posterior mean of Beta distributoin 
#' * Posterior SD - posterior standard deviation of Beta distributoin 
#' * LL - lower limit of the credible interval
#' * UL - upper limit of the credible interval
#'
#'
#' @references
#' \insertRef{Gelman2004}{statpsych}                            
#'
#'
#' @examples
#' ci.bayes.prop(.05, .4, .1, 12, 100)
#'
#' # Should return:
#' # Posterior mean Posterior SD       LL        UL
#' #           0.15   0.03273268  0.09218 0.2188484
#'
#'
#' @importFrom stats qnorm
#' @importFrom stats qbeta
#' @export
ci.bayes.prop <- function(alpha, prior.mean, prior.sd, f, n) {
 if (prior.sd^2 >= prior.mean*(1 - prior.mean)) {stop("prior SD is too large")}
 if (f > n) {stop("f cannot be greater than n")}
 zcrit <- qnorm(1 - alpha/2)
 a <- ((1 - prior.mean)/prior.sd^2 - 1/prior.mean)*prior.mean^2
 b <- a*(1/prior.mean - 1)
 post.mean <- (f + a)/(n + a + b)
 post.sd <- sqrt((f + a)*(n - f + b)/((n + a + b)^2*(a + b + n - 1)))
 post.a <- a + f
 post.b <- b + n - f
 ll <- qbeta(alpha/2, post.a, post.b)
 ul <- qbeta(1 - alpha/2, post.a, post.b)
 out <- t(c(post.mean, post.sd, ll, ul))
 colnames(out) <- c("Posterior mean", "Posterior SD", "LL", "UL")
 rownames(out) <- ""
 return(out)
}


#  ci.pv =====================================================================
#' Confidence intervals for positive and negative predictive values with 
#' retrospective sampling
#'
#'                                 
#' @description
#' Computes adjusted Wald confidence intervals for positive and negative
#' predictive values (PPV and NPV) of a diagnostic test with retrospective 
#' sampling where the population prevalence rate is assumed to be known. With 
#' retrospective sampling, one random sample is obtained from a subpopulation
#' that is known to have a "positive" outcome, a second random sample is
#' obtained from a subpopulation that is known to have a "negative" outcome,
#' and then the diagnostic test (scored "pass" or "fail") is given in each 
#' sample. PPV and NPV can be expressed as a function of proportion ratios 
#' and the known population prevalence rate (the population proportion who 
#' would "pass"). The confidence intervals for PPV and NPV are based on the 
#' Price-Bonett adjusted Wald confidence interval for a proportion ratio.
#'
#'
#' @param   alpha   alpha level for 1-alpha confidence
#' @param   f1      number of participants with a positive outcome who pass the test
#' @param   f2      number of participants with a negative outcome who fail the test
#' @param   n1      sample size for the positive outcome group
#' @param   n2      sample size for the negative outcome group
#' @param   prev    known population proportion with a positive outcome
#'
#'
#' @return
#' Returns a 2-row matrix. The columns are:
#' * Estimate - adjusted estimate of the predictive value
#' * LL - lower limit of the adjusted Wald confidence interval
#' * UL - upper limit of the adjusted Wald confidence interval
#'
#'
#' @references
#' \insertRef{Price2008}{statpsych}
#'
#'
#' @examples
#' ci.pv(.05, 89, 5, 100, 100, .16)
#'
#' # Should return:
#' #        Estimate        LL        UL
#' # PPV:  0.7640449 0.5838940 0.8819671
#' # NPV:  0.9779978 0.9623406 0.9872318
#'
#'
#' @importFrom stats qnorm
#' @export
ci.pv <- function(alpha, f1, f2, n1, n2, prev) {
 z <- qnorm(1 - alpha/2)
 k1 <- (1 - prev)/prev
 p1 <- (f1 + 1/4)/(n1 + 7/4)
 p2 <- (f2 + 1/4)/(n2 + 7/4)
 v1 <- 1/(f1 + 1/4 + (f1 + 1/4)^2/(n1 - f1 + 3/2))
 v2 <- 1/(f2 + 1/4 + (f2 + 1/4)^2/(n2 - f2 + 3/2))
 se <- sqrt(v1 + v2)
 LL0 <- exp(log(p2/p1) - z*se)
 UL0 <- exp(log(p2/p1) + z*se)
 ppv <- 1/(1 + (p2/p1)*k1)
 LL1 <- 1/(1 + UL0*k1)
 UL1 <- 1/(1 + LL0*k1)
 out1 <- t(c(ppv, LL1, UL1))
 k2 <- prev/(1 - prev)
 f3 <- n1 - f1
 f4 <- n2 - f2
 p3 <- (f3 + 1/4)/(n1 + 7/4)
 p4 <- (f4 + 1/4)/(n2 + 7/4)
 v1 <- 1/(f3 + 1/4 + (f3 + 1/4)^2/(n1 - f3 + 3/2))
 v2 <- 1/(f2 + 1/4 + (f4 + 1/4)^2/(n2 - f4 + 3/2))
 se <- sqrt(v1 + v2)
 LL0 <- exp(log(p3/p4) - z*se)
 UL0 <- exp(log(p3/p4) + z*se)
 npv <- 1/(1 + (p3/p4)*k2)
 LL2 <- 1/(1 + UL0*k2)
 UL2 <- 1/(1 + LL0*k2)
 out2 <- t(c(npv, LL2, UL2))
 out <- rbind(out1, out2)
 colnames(out) <- c("Estimate", "LL", "UL")
 rownames(out) <- c("PPV: ", "NPV: ")
 return(out)
}


#  ci.poisson =============================================================== 
#' Confidence interval for a Poisson rate
#'
#'                        
#' @description
#' Computes a confidence interval for a population Poisson rate. This function
#' requires the number of occurances (f) of a specific event that were 
#' observed over a specific period of time (t).
#'
#'
#' @param  alpha  alpha value for 1-alpha confidence 
#' @param  f      number of event occurances
#' @param  t      time period 
#'
#'
#' @details
#' The time period (t) does not need to be an integer and can be expressed in 
#' any unit of time such as seconds, hours, or months. The occurances are
#' assumed to be independent of one another and the unknown occurance rate is
#' assumed to be constant over time. 
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimated Poisson rate
#' * SE - recovered standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @references
#' \insertRef{Hahn1991}{statpsych}
#'
#'
#' @examples
#' ci.poisson(.05, 23, 5.25)
#'
#' # Should return:
#' # Estimate        SE       LL      UL
#' # 4.380952 0.9684952 2.777148 6.57358

#'  
#' 
#' @importFrom stats qchisq
#' @importFrom stats qnorm
#' @export
ci.poisson <- function(alpha, f, t) {
 z <- qnorm(1 - alpha/2)
 est <- f/t
 ll <- qchisq(alpha/2, 2*f)/(2*t)
 ul <- qchisq(1 - alpha/2, 2*f + 2)/(2*t)
 se <- (ul - ll)/(2*z)
 out <- t(c(est, se, ll, ul))
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 rownames(out) <- ""
 return(out)
}


#  ci.ratio.poisson2 ========================================================== 
#' Confidence interval for a ratio of Poisson rates in a 2-group design
#'
#'                        
#' @description
#' Computes a confidence interval for a ratio of population Poisson rates in a
#' 2-group design. The confidence interval is based on the binomial method 
#' with an Agresti-Coull confidence interval. This function requires the number 
#' of occurances of a specific event (f) that were observed over a specific
#' period of time (t) within each group.
#'
#'
#' @param  alpha  alpha value for 1-alpha confidence 
#' @param  f1     number of event occurances for group 1
#' @param  f2     number of event occurances for group 2
#' @param  t1     time period for group 1
#' @param  t2     time period for group 2
#'
#'
#' @details
#' The time periods do not need to be integers and can be expressed in any unit
#' of time such as seconds, hours, or months. The occurances are assumed to be
#' independent of one another and the unknown occurance rate is assumed to be
#' constant over time within each group condition.
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimated ratio of Poisson rates
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @references
#' \insertRef{Price2000}{statpsych}
#'
#'
#' @examples
#' ci.ratio.poisson2(.05, 19, 5, 30, 40.5)
#'
#' # Should return:
#' # Estimate       LL       UL
#' #     5.13 1.939576 13.71481
#'  
#' 
#' @importFrom stats qnorm
#' @export
ci.ratio.poisson2 <- function(alpha, f1, f2, t1, t2) {
 z <- qnorm(1 - alpha/2)
 est <- (f1/t1)/(f2/t2)
 n <- f1 + f2
 p <- (f1 + 2)/(n + 4)
 se <- sqrt(p*(1 - p)/(n + 4))
 ll0 <- (p - z*se)
 ul0 <- (p + z*se)
 ll <- (t2/t1)*ll0/(1 - ll0)
 ul <- (t2/t1)*ul0/(1 - ul0)
 out <- t(c(est, ll, ul))
 colnames(out) <- c("Estimate", "LL", "UL")
 rownames(out) <- ""
 return(out)
}


#  pi.prop =================================================================== 
#' Prediction interval for an estimated proportion 
#'
#'                        
#' @description
#' Computes approximate prediction interval for the estimated proportion 
#' in a future study with a planned sample size of n2. The prediction interval
#' uses a proportion estimate from a prior study that used a sample size of n1.
#'
#'
#' @param  alpha  alpha value for 1-alpha confidence 
#' @param  prop   estimated proportion from prior study
#' @param  n1     sample size used to estimate proportion in prior study 
#' @param  n2     planned sample size of future study
#'
#'
#' @return 
#' Returns a prediction interval for an estimated proportion in a future 
#' study
#'
#'
#' @examples
#' pi.prop(.1, .225, 80, 120)
#'
#' # Should return:
#' #         LL       UL
#' #  0.1390955 0.337095
#'  
#' 
#' @importFrom stats qnorm
#' @export
pi.prop <- function(alpha, prop, n1, n2) {
 z <- qnorm(1 - alpha/2)
 p <- (n1*prop + 2)/(n1 + 4)
 ll <- p - z*sqrt(p*(1 - p)/(n1 + 4) + p*(1 - p)/(n2 + 4))
 ul <- p + z*sqrt(p*(1 - p)/(n1 + 4) + p*(1 - p)/(n2 + 4))
 out <- t(c(ll, ul))
 colnames(out) <- c("LL", "UL")
 rownames(out) <- ""
 return(out)
}


# ======================== Hypothesis Tests ==================================
# test.prop =================================================================
#' Hypothesis test for a single proportion 
#'
#'
#' @description
#' Computes a continuity-corrected z-test for a single proportion in a 
#' 1-group design. A confidence interval for a population proportion 
#' is a recommended supplement to the z-test.
#'
#'
#' @param   f    number of participants who have the attribute
#' @param   n    sample size
#' @param   h    null hypothesis value of proportion
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - ML estimate of proportion 
#' * z - z test statistic
#' * p - two-sided p-value
#'
#'
#' @references
#' \insertRef{Snedecor1980}{statpsych}
#'
#'
#' @examples
#' test.prop(9, 20, .2)
#'
#' # Should return:
#' # Estimate        z          p
#' #     0.45 2.515576 0.01188379
#'
#'
#' @importFrom stats pnorm
#' @export
test.prop <- function(f, n, h) {
 if (f > n) {stop("f cannot be greater than n")}
 p <- f/n
 se <- sqrt(h*(1 - h)/n)
 z <- (abs(p - h) - 1/(2*n))/se
 pval <- 2*(1 - pnorm(abs(z)))
 out <- t(c(p, z, pval))
 colnames(out) <- c("Estimate", "z", "p")
 rownames(out) <- ""
 return(out)
}


# test.prop2 =================================================================
#' Hypothesis test for a 2-group proportion difference
#'
#'
#' @description
#' Computes a continuity-corrected z-test for a difference of proportions in  
#' a 2-group design. A confidence interval for a difference in population
#' proportions (see ci.prop2) is a recommended supplement to the z-test.
#'
#'
#' @param  f1      number of group 1 participants who have the attribute
#' @param  f2      number of group 2 participants who have the attribute
#' @param  n1      sample size for group 1
#' @param  n2      sample size for group 2
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - ML estimate of proportion difference
#' * z - z test statistic
#' * p - two-sided p-value
#'
#'
#' @references
#' \insertRef{Snedecor1980}{statpsych}
#'
#'
#' @examples
#' test.prop2(11, 26, 50, 50)
#'
#' # Should return:
#' # Estimate        z           p
#' #     -0.3 2.899726 0.003734895
#'
#'
#' @importFrom stats pnorm
#' @export
test.prop2 <- function(f1, f2, n1, n2) {
 if (f1 > n1) {stop("f cannot be greater than n")}
 if (f2 > n2) {stop("f cannot be greater than n")}
 p1 <- f1/n1
 p2 <- f2/n2
 diff <- p1 - p2
 p0 <- (f1 + f2)/(n1 + n2)
 se <- sqrt(p0*(1 - p0)/n1 + p0*(1 - p0)/n2)
 z <- (abs(p1 - p2) - 1/(2*n1) - 1/(2*n2))/se
 pval <- 2*(1 - pnorm(abs(z)))
 out <- t(c(diff, z, pval))
 colnames(out) <- c("Estimate", "z", "p")
 rownames(out) <- ""
 return(out)
}


#  test.prop.bs ==============================================================
#' Hypothesis test of equal proportions in a between-subjects design
#'
#'
#' @description
#' Computes a Pearson chi-square test for equal population proportions for a 
#' dichotomous response variable in a one-factor between-subjects design.
#'
#'
#' @param  f    vector of frequency counts of participants who have the attribute
#' @param  n    vector of sample sizes
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Chi-square - chi-square test statistic
#' * df - degrees of freedom 
#' * p - p-value
#'
#'
#' @references
#' \insertRef{Fleiss2003}{statpsych}
#'
#'
#' @examples
#' f <- c(35, 30, 15)
#' n <- c(50, 50, 50)
#' test.prop.bs (f, n)
#'
#' # Should return:
#' # Chi-square df            p
#' #   17.41071  2 0.0001656958
#'
#'
#' @importFrom stats pchisq
#' @export
test.prop.bs <- function(f, n) {
 s <- sum(as.integer(as.logical(n < f)))
 if (s > 0) {stop("f cannot be greater than n")}
 p0 <- sum(f)/sum(n)
 df <- length(f) - 1
 a <- 1/(p0*(1 - p0))
 p <- f/n
 r <- (p - p0)^2
 chi <- a*sum(n*r)
 pval <- 1 - pchisq(chi, df)
 out <- t(c(chi, df, pval))
 colnames(out) <- c("Chi-square", "df", "p")
 rownames(out) <- ""
 return(out)
}


#  test.prop.ps ============================================================== 
#' Hypothesis test for a paired-samples proportion difference
#'
#'
#' @description
#' Computes a continuity-corrected McNemar test for equality of proportions 
#' in a paired-samples design. This function requires the frequency counts
#' from a 2 x 2 contingency table for two paired dichotomous measurements.
#' A confidence interval for a difference in population proportions (see
#' ci.prop.ps) is a recommended supplement to the McNemar test.
#'
#'
#' @param   f00    number participants with y = 0 and x = 0
#' @param   f01    number participants with y = 0 and x = 1
#' @param   f10    number participants with y = 1 and x = 0
#' @param   f11    number participants with y = 1 and x = 1
#'
#'
#' @references
#' \insertRef{Snedecor1980}{statpsych}
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - ML estimate of proportion difference
#' * z - z test statistic
#' * p - two-sided p-value
#'
#'
#' @examples
#' test.prop.ps(156, 96, 68, 80)
#'
#' # Should return:
#' # Estimate        z          p
#' #     0.07 2.108346 0.03500109
#'
#'
#' @importFrom stats pnorm
#' @export 
test.prop.ps <- function(f00, f01, f10, f11) {
 n <- f00 + f01 + f10 + f11
 p1 <- (f01 + f11)/n
 p2 <- (f10 + f11)/n
 p01 <- f01/n
 p10 <- f10/n
 diff <- p1 - p2
 se <- sqrt((p10 + p01)/n)
 z <- (abs(p10 - p01) - 1/n)/se
 pval <- 2*(1 - pnorm(abs(z)))
 out <- t(c(diff, z, pval))
 colnames(out) <- c("Estimate", "z", "p")
 rownames(out) <- ""
 return(out)
}


#  test.mono.prop.bs ============================================================
#' Test of monotonic trend in proportions for an ordered between-subjects
#' factor
#'
#'
#' @description
#' Computes simultaneous confidence intervals for all adjacent pairwise
#' comparisons of population proportions using group frequency counts and
#' samples sizes as input. If one or more lower limits are greater than
#' 0 and no upper limit is less than 0, then conclude that the population
#' proportions are monotoic decreasing. If one or more upper limits are 
#' less than 0 and no lower limits are greater than 0, then conclude that
#' the population proportions are monotoic increasing. Reject the hypothesis
#' of a monotonic trend if any lower limit is greater than 0 and any upper 
#' limit is less than 0. 
#'
#'
#' @param  alpha   alpha level for simultaneous 1-alpha confidence
#' @param  f       vector of frequency counts of participants who have the attribute
#' @param  n       vector of sample sizes
#'
#'
#' @return 
#' Returns a matrix with the number of rows equal to the number
#' of adjacent pairwise comparisons. The columns are:
#' * Estimate - estimated proportion difference
#' * SE - standard error
#' * LL - one-sided lower limit of the confidence interval
#' * UL - one-sided upper limit of the confidence interval
#'
#'
#' @examples
#' f <- c(67, 49, 30, 10)
#' n <- c(100, 100, 100, 100)
#' test.mono.prop.bs(.05, f, n)
#'
#' # Should return:
#' #      Estimate         SE         LL        UL
#' # 1 2 0.1764706 0.06803446 0.01359747 0.3393437
#' # 2 3 0.1862745 0.06726135 0.02525219 0.3472968
#' # 3 4 0.1960784 0.05493010 0.06457688 0.3275800
#'
#'
#' @importFrom stats qnorm
#' @export
test.mono.prop.bs <-function(alpha, f, n) {
 s <- sum(as.integer(as.logical(n < f)))
 if (s > 0) {stop("f cannot be greater than n")}
 a <- length(f)
 p.adj <- (f + 1)/(n + 2)
 v <- p.adj*(1 - p.adj)/(n + 2)
 p1 <- p.adj[1: a - 1]
 p2 <- p.adj[2: a]
 Estimate <- p1 - p2
 v1 <- v[1: a - 1]
 v2 <- v[2: a]
 n1 <- n[1: a - 1]
 n2 <- n[2: a]
 SE <- sqrt(v1 + v2)
 zcrit <- qnorm(1 - alpha/(2*(a - 1)))
 LL <- Estimate - zcrit*SE
 UL <- Estimate + zcrit*SE
 pair <- cbind(seq(1, a - 1), seq(2, a))
 out <- cbind(pair, Estimate, SE, LL, UL)
 rownames(out) <- rep("", a - 1)
 return(out)
}


# ================= Sample Size for Desired Precision =======================
#  size.ci.prop ============================================================
#' Sample size for a single proportion confidence interval  
#'
#'
#' @description
#' Computes the sample size required to estimate a single proportion with 
#' desired confidence interval precision. Set the proportion planning value
#' to .5 for a conservatively large sample size.
#'
#'
#' @param  alpha  alpha level for 1-alpha confidence 
#' @param  p      planning value of proportion
#' @param  w      desired confidence interval width
#'
#'
#' @return
#' Returns the required sample size
#'
#'
#' @examples
#' size.ci.prop(.05, .4, .2)
#'
#' # Should return:
#' # Sample size
#' #          93
#'
#'
#' @importFrom stats qnorm
#' @export                 
size.ci.prop <- function(alpha, p, w) {
 if (p > .9999 || p < .0001) {stop("proportion must be between .0001 and .9999")}
 z <- qnorm(1 - alpha/2)
 n <- ceiling(4*p*(1 - p)*(z/w)^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 rownames(out) <- ""
 return(out)
}


#  size.ci.prop2 =============================================================
#' Sample size for a 2-group proportion difference confidence interval  
#'
#'
#' @description
#' Computes the sample size in each group (assuming equal sample sizes) required
#' to estimate a difference of proportions with desired confidence interval 
#' precision in a 2-group design. Set the proportion planning values to .5 for
#' a conservatively large sample size.
#'
#'
#' @param  alpha  alpha level for 1-alpha confidence 
#' @param  p1     planning value of proportion for group 1
#' @param  p2     planning value of proportion for group 2
#' @param  w      desired confidence interval width
#'
#'
#' @return
#' Returns the required sample size for each group
#'
#'
#' @examples
#' size.ci.prop2(.05, .4, .2, .15)
#'
#' # Should return:
#' # Sample size per group
#' #                   274
#'
#'
#' @importFrom stats qnorm
#' @export                 
size.ci.prop2 <- function(alpha, p1, p2, w) {
 if (p1 > .9999 || p1 < .0001) {stop("p1 must be between .0001 and .9999")}
 if (p2 > .9999 || p2 < .0001) {stop("p2 must be between .0001 and .9999")}
 z <- qnorm(1 - alpha/2)
 n <- ceiling(4*(p1*(1 - p1) + p2*(1 - p2))*(z/w)^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size per group"
 rownames(out) <- ""
 return(out)
}


#  size.ci.ratio.prop2 =======================================================
#' Sample size for a 2-group proportion ratio confidence interval  
#'
#'
#' @description
#' Computes the sample size in each group (assuming equal sample sizes) required
#' to estimate a ratio of proportions with desired confidence interval precision 
#' in a 2-group design.
#'
#'
#' @param  alpha  alpha level for 1-alpha confidence 
#' @param  p1     planning value of proportion for group 1
#' @param  p2     planning value of proportion for group 2
#' @param  r      desired upper to lower confidence interval endpoint ratio
#'
#'
#' @return
#' Returns the required sample size for each group
#'
#'
#' @examples
#' size.ci.ratio.prop2(.05, .2, .1, 2)
#'
#' # Should return:
#' # Sample size per group
#' #                   416
#'
#'
#' @importFrom stats qnorm
#' @export   
size.ci.ratio.prop2 <- function(alpha, p1, p2, r) {
 if (p1 > .9999 || p1 < .0001) {stop("p1 must be between .0001 and .9999")}
 if (p2 > .9999 || p2 < .0001) {stop("p2 must be between .0001 and .9999")}
 z <- qnorm(1 - alpha/2)
 n <- ceiling(4*((1 - p1)/p1 + (1 - p2)/p2)*(z/log(r))^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size per group"
 rownames(out) <- ""
 return(out)
}


#  size.ci.lc.prop.bs =========================================================
#' Sample size for a between-subjects proportion linear contrast confidence 
#' interval  
#'
#'
#' @description
#' Computes the sample size in each group (assuming equal sample sizes) required 
#' to estimate a linear contrast of proportions with desired confidence interval 
#' precision in a between-subjects design. Set the proportion planning values to
#' .5 for a conservatively large sample size.
#'
#'
#' @param  alpha  alpha level for 1-alpha confidence 
#' @param  p      vector of proportion planning values
#' @param  w      desired confidence interval width
#' @param  v      vector of between-subjects contrast coefficients
#'
#'
#' @return
#' Returns the required sample size for each group
#'
#'
#' @examples
#' p <- c(.25, .30, .50, .50)
#' v <- c(.5, .5, -.5, -.5)
#' size.ci.lc.prop.bs(.05, p, .2, v)
#'
#' # Should return:
#' # Sample size per group
#' #                    87
#'
#'
#' @importFrom stats qnorm
#' @export
size.ci.lc.prop.bs <- function(alpha, p, w, v) {
 z <- qnorm(1 - alpha/2)
 n <- ceiling((4*t(v)%*%diag(p*(1 - p))%*%v)*(z/w)^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size per group"
 rownames(out) <- ""
 return(out)
}


#  size.ci.prop.ps ===========================================================
#' Sample size for a paired-sample proportion difference confidence interval  
#'
#'
#' @description
#' Computes the sample size required to estimate a proportion difference with 
#' desired confidence interval precision in a paired-samples design. Set the 
#' proportion planning values to .5 for a conservatively large sample size.
#' Set the phi correlation planning value to the smallest value within a
#' plausible range for a conservatively large sample size.
#'
#'
#' @param  alpha  alpha level for 1-alpha confidence 
#' @param  p1     planning value of proportion for measurement 1
#' @param  p2     planning value of proportion for measurement 2
#' @param  phi    planning value of phi correlation
#' @param  w      desired confidence interval width
#'
#'
#' @return
#' Returns the required sample size
#'
#'
#' @examples
#' size.ci.prop.ps(.05, .2, .3, .8, .1)
#'
#' # Should return:
#' # Sample size
#' #         118
#'
#'
#' @importFrom stats qnorm
#' @export  
size.ci.prop.ps <- function(alpha, p1, p2, phi, w) {
 if (phi > .999 || phi < -.999) {stop("phi must be between -.999 and .999")}
 if (p1 > .9999 || p1 < .0001) {stop("p1 must be between .0001 and .9999")}
 if (p2 > .9999 || p2 < .0001) {stop("p2 must be between .0001 and .9999")}
 z <- qnorm(1 - alpha/2)
 cov <- phi*sqrt(p1*p2*(1 - p1)*(1 - p2))
 n <- ceiling((4*(p1*(1 - p1) + p2*(1 - p2) - 2*cov))*(z/w)^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 rownames(out) <- ""
 return(out)
}


#  size.ci.ratio.prop.ps ======================================================
#' Sample size for a paired-samples proportion ratio confidence interval  
#'
#'
#' @description
#' Computes the sample size required to estimate a ratio of proportions 
#' with desired confidence interval precision in a paired-samples design.
#' Set the phi correlation planning value to the smallest value within a 
#' plausible range for a conservatively large sample size. 
#'
#'
#' @param  alpha  alpha level for 1-alpha confidence 
#' @param  p1     planning value of proportion for measurement 1
#' @param  p2     planning value of proportion for measurement 2
#' @param  phi    planning value of phi correlation
#' @param  r      desired upper to lower confidence interval endpoint ratio
#'
#'
#' @return
#' Returns the required sample size
#'
#'
#' @examples
#' size.ci.ratio.prop.ps(.05, .4, .2, .7, 2)
#'
#' # Should return:
#' # Sample size
#' #          67
#'
#'
#' @importFrom stats qnorm
#' @export  
size.ci.ratio.prop.ps <- function(alpha, p1, p2, phi, r) {
 if (phi > .999 || phi < -.999) {stop("phi must be between -.999 and .999")}
 if (p1 > .9999 || p1 < .0001) {stop("p1 must be between .0001 and .9999")}
 if (p2 > .9999 || p2 < .0001) {stop("p2 must be between .0001 and .9999")}
 z <- qnorm(1 - alpha/2)
 cov <- phi*sqrt((1 - p1)*(1 - p2)/(p1*p2))
 n <- ceiling(4*((1 - p1)/p1 + (1 - p2)/p2 - 2*cov)*(z/log(r))^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 rownames(out) <- ""
 return(out)
}


#  size.ci.agree =============================================================
#' Sample size for a G-index confidence interval
#'
#'
#' @description
#' Computes the sample size required to estimate a G-index of agreement for 
#' two dichotomous ratings with desired confidence interval precision.
#' Set the G-index planning value to the smallest value within a plausible 
#' range for a conservatively large sample size.
#'
#'
#' @param  alpha  alpha level for 1-alpha confidence 
#' @param  G      planning value of G-index
#' @param  w      desired confidence interval width
#'
#'
#' @return
#' Returns the required sample size
#'
#'
#' @examples
#' size.ci.agree(.05, .8, .2)
#'
#' # Should return:
#' # Sample size
#' #         139
#'
#'
#' @importFrom stats qnorm
#' @export
size.ci.agree <- function(alpha, G, w) {
 if (G > .999 || G < .001) {stop("G must be between .001 and .999")}
 z <- qnorm(1 - alpha/2)
 n <- ceiling(4*(1 - G^2)*(z/w)^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 rownames(out) <- ""
 return(out)
}


# ===================== Sample Size for Desired Power ========================
#  size.test.prop =========================================================== 
#' Sample size for a test of a single proportion 
#'
#' @description
#' Computes the sample size required to test a single population proportion 
#' with desired power using a correction for continuity in a 1-group design. 
#'
#'
#' @param  alpha  alpha level for hypothesis test 
#' @param  pow    desired power
#' @param  p      planning value of proportion 
#' @param  h      null hypothesis value of proportion
#'
#'
#' @references
#' \insertRef{Fleiss2003}{statpsych}
#'
#'
#' @return
#' Returns the required sample size
#'
#'
#' @examples
#' size.test.prop(.05, .9, .5, .3)
#'
#' # Should return:
#' # Sample size
#' #          67
#'
#'
#' @importFrom stats qnorm
#' @export
size.test.prop <- function(alpha, pow, p, h) {
 if (p > .9999 || p < .0001) {stop("proportion must be between .0001 and .9999")}
 if (h > .9999 || h < .0001) {stop("null hypothesis value must be between .0001 and .9999")}
 za <- qnorm(1 - alpha/2)
 zb <- qnorm(pow)
 n0 <- ceiling((za*sqrt(h*(1 - h)) + zb*sqrt(p*(1 - p)))^2/((p - h)^2))
 n <- n0 + 1/abs(p - h)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 rownames(out) <- ""
 return(out)
}


#  size.test.prop2 ===========================================================
#' Sample size for a test of a 2-group proportion difference
#'
#'
#' @description
#' Computes the sample size in each group required to test a difference in 
#' population proportions with desired power and a continuity correction in a
#' 2-group design. This function requires planning values for both proportions. 
#' Set each proportion planning value to .5, or a value closest to .5 within
#' a plausible range, for a conservatively large sample size requirement. This
#' function does not use the typical sample size approach where the effect 
#' size is assumed to equal the difference in proportion planning values. This
#' function does not require the planning value for the proportion difference 
#' (effect size) to equal the difference of the two proportion planning values;
#' for example, the planning value of the proportion difference could be set 
#' equal to a minimally interesting effect size.
#'
#'
#' @param  alpha  alpha level for hypothesis test 
#' @param  pow    desired power
#' @param  p1     planning value of proportion for group 1
#' @param  p2     planning value of proportion for group 2
#' @param  es     planning value of proportion difference (effect size)
#'
#'
#' @return
#' Returns the required sample size for eachr group
#'
#'
#' @examples
#' size.test.prop2(.05, .8, .5, .5, .2)
#'
#' # Should return:
#' # Sample size per group
#' #                   109
#'
#' size.test.prop2(.05, .8, .3, .1, .2)
#' # Should return:
#' # Sample size per group
#' #                    71
#'
#'
#' @importFrom stats qnorm
#' @export
size.test.prop2 <- function(alpha, pow, p1, p2, es) {
 if (p1 > .9999 || p1 < .0001) {stop("p1 must be between .0001 and .9999")}
 if (p2 > .9999 || p2 < .0001) {stop("p2 must be between .0001 and .9999")}
 if (es < .001) {stop("effect size must be greater than .001")}
 za <- qnorm(1 - alpha/2)
 zb <- qnorm(pow)
 p0 <- (p1 + p2)/2
 se1 <- sqrt(2*(p0*(1 - p0)))
 se2 <- sqrt(p1*(1 - p1) + p2*(1 - p2))
 n0 <- ceiling((za*se2 + zb*se1)^2/es^2)
 n <- n0 + 2/abs(es)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size per group"
 rownames(out) <- ""
 return(out)
}


#  size.test.lc.prop.bs ======================================================= 
#' Sample size for a test of between-subjects proportion linear contrast
#'
#'
#' @description
#' Computes the sample size in each group (assuming equal sample sizes) required
#' to test a linear contrast of population proportions with desired power in a 
#' between-subjects design. This function requires planning values for all 
#' proportions. Set the proportion planning values to .5 for a conservatively 
#' large sample size. The planning value for the effect size (linear contrast of 
#' proportions) could be set equal to the linear contrast of proportion planning 
#' values or it could be set equal to a minimally interesting effect size.
#'
#'
#' @param  alpha  alpha level for hypothesis test 
#' @param  pow    desired power
#' @param  p      vector of proportion planning values
#' @param  es     planning value of proportion linear contrast
#' @param  v      vector of between-subjects contrast coefficients
#'
#'
#' @return
#' Returns the required sample size for each group
#'
#'
#' @examples
#' p <- c(.25, .30, .50, .50)
#' v <- c(.5, .5, -.5, -.5)
#' size.test.lc.prop.bs(.05, .9, p, .15, v)
#'
#' # Should return:
#' # Sample size per group
#' #                   105
#'
#'
#' @importFrom stats qnorm
#' @export
size.test.lc.prop.bs <- function(alpha, pow, p, es, v) {
 za <- qnorm(1 - alpha/2)
 zb <- qnorm(pow)
 n <- ceiling((t(v)%*%diag(p*(1 - p))%*%v)*(za + zb)^2/es^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size per group"
 rownames(out) <- ""
 return(out)
}


#  size.equiv.prop2 ==========================================================
#' Sample size for a 2-group proportion equivalence test
#'
#'
#' @description
#' Computes the sample size in each group (assuming equal sample sizes) required
#' to perform an equivalence test for the difference in population proportions 
#' with desired power in a 2-group design. The value of h specifies a range of
#' practical equivalence, -h to h, for the difference in population proportions. 
#' The absolute difference in the proportion planning values must be less than h.  
#' Equivalence tests often require a very large sample size. Equivalence tests 
#' usually use 2 x alpha rather than alpha (e.g., use alpha = .10 rather than
#' alpha = .05).
#'
#'
#' @param  alpha  alpha level for hypothesis test 
#' @param  pow    desired power
#' @param  p1     planning value of proportion for group 1
#' @param  p2     planning value of proportion for group 2
#' @param  h      upper limit for range of practical equivalence
#'
#'
#' @return
#' Returns the required sample size for each group
#'
#'
#' @examples
#' size.equiv.prop2(.1, .8, .30, .35, .15)
#'
#' # Should return:
#' # Sample size per group
#' #                   288
#'
#'
#' @importFrom stats qnorm
#' @export
size.equiv.prop2 <- function(alpha, pow, p1, p2, h) {
 if (h <= abs(p1 - p2)) {stop("|p1 - p2| must be less than h")}
 if (p1 > .9999 || p1 < .0001) {stop("p1 must be between .0001 and .9999")}
 if (p2 > .9999 || p2 < .0001) {stop("p2 must be between .0001 and .9999")}
 za <- qnorm(1 - alpha)
 zb <- qnorm(1 - (1 - pow)/2)
 es <- p1 - p2
 n <- ceiling((p1*(1 - p1) + p2*(1 - p2))*(za + zb)^2/(h - abs(es))^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size per group"
 rownames(out) <- ""
 return(out)
}


#  size.supinf.prop2 ==========================================================
#' Sample size for a 2-group superiority or inferiority test of proportions
#'
#'
#' @description
#' Computes the sample size in each group (assuming equal sample sizes) required
#' to perform a superiority or inferiority test for the difference in population 
#' proportions with desired power in a 2-group design. For a superiority test, 
#' specify the upper limit (h) for the range of practical equivalence and specify
#' values of p1 and p2 such that p1 - p2 > h. For an inferiority test, specify the 
#' lower limit (-h) for the range of practical equivalence and specify values
#' of p1 and p2 such that p1 - p2 > -h. This function sets the effect size equal
#' to p1 - p2.
#'
#'
#' @param  alpha  alpha level for hypothesis test 
#' @param  pow    desired power
#' @param  p1     planning value of proportion for group 1
#' @param  p2     planning value of proportion for group 2
#' @param  h      lower or upper limit for range of practical equivalence
#'
#'
#' @return
#' Returns the required sample size for each group
#'
#'
#' @examples
#' size.supinf.prop2(.05, .9, .35, .20, .05)
#'
#' # Should return:
#' # Sample size per group
#' #                   408
#'
#'
#' @importFrom stats qnorm
#' @export
size.supinf.prop2 <- function(alpha, pow, p1, p2, h) {
 if (p1 > .9999 || p1 < .0001) {stop("p1 must be between .0001 and .9999")}
 if (p2 > .9999 || p2 < .0001) {stop("p2 must be between .0001 and .9999")}
 za <- qnorm(1 - alpha/2)
 zb <- qnorm(pow)
 es <- p1 - p2
 n <- ceiling((p1*(1 - p1) + p2*(1 - p2))*(za + zb)^2/(es - h)^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size per group"
 rownames(out) <- ""
 return(out)
}


#  size.test.prop.ps =========================================================
#' Sample size for a test of a paired-samples proportion difference
#'
#'
#' @description
#' Computes the sample size required to test a difference in population 
#' proportions with desired power in a paired-samples design. This function 
#' requires planning values for both proportions and a phi coefficient that 
#' describes the correlation between the two dichotomous measurements. The 
#' proportion planning values can be set to .5 for a conservatively large 
#' sample size. The planning value for the effect size (proportion difference)
#' could be set equal to the difference of the two proportion planning values 
#' or it could be set equal to a minimally interesting effect size. Set the
#' phi correlation planning value to the smallest value within a plausible range
#' for a conservatively large sample size.
#'
#'
#' @param  alpha  alpha level for hypothesis test 
#' @param  pow    desired power
#' @param  p1     planning value of proportion for measurement 1
#' @param  p2     planning value of proportion for measurement 2
#' @param  phi    planning value of phi correlation
#' @param  es     planning value of proportion difference
#'
#'
#' @return
#' Returns the required sample size
#'
#'
#' @examples
#' size.test.prop.ps(.05, .80, .4, .3, .5, .1)
#'
#' # Should return:
#' # Sample size
#' #         177
#'
#'
#' @importFrom stats qnorm
#' @export
size.test.prop.ps <- function(alpha, pow, p1, p2, phi, es) {
 if (phi > .999 || phi < -.999) {stop("phi must be between -.999 and .999")}
 za <- qnorm(1 - alpha/2)
 zb <- qnorm(pow)
 cov <- phi*sqrt(p1*p2*(1 - p1)*(1 - p2))
 n <- ceiling(((p1*(1 - p1) + p2*(1 - p2) - 2*cov))*((za + zb)/es)^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 rownames(out) <- ""
 return(out)
}


#  size.equiv.prop.ps ========================================================
#' Sample size for a paired-samples proportion equivalence test
#'
#'
#' @description
#' Computes the sample size required to perform an equivalence test for the
#' difference in population proportions with desired power in a paired-samples
#' design. The value of h specifies a range of practical equivalence, -h to h, 
#' for the difference in population proportions. The absolute difference in 
#' the proportion planning values must be less than h.  Equivalence tests often
#' require a very large sample size. Equivalence tests usually use 2 x alpha 
#' rather than alpha (e.g., use alpha = .10 rather alpha = .05). This function
#' sets the effect size equal to the difference in proportion planning values. 
#' Set the phi correlation planning value to the smallest value within a 
#' plausible range for a conservatively large sample size.
#'
#'
#' @param  alpha  alpha level for hypothesis test 
#' @param  pow    desired power
#' @param  p1     planning value of proportion for measurement 1
#' @param  p2     planning value of proportion for measurement 2
#' @param  phi    planning value of phi correlation
#' @param  h      upper limit for range of practical equivalence
#'
#'
#' @return
#' Returns the required sample size
#'
#'
#' @examples
#' size.equiv.prop.ps(.1, .8, .30, .35, .40, .15)
#'
#' # Should return:
#' # Sample size
#' #         173
#'
#'
#' @importFrom stats qnorm
#' @export
size.equiv.prop.ps <- function(alpha, pow, p1, p2, phi, h) {
 if (phi > .999 || phi < -.999) {stop("phi must be between -.999 and .999")}
 if (h <= abs(p1 - p2)) {stop("|p1 - p2| must be less than h")}
 if (p1 > .9999 || p1 < .0001) {stop("p1 must be between .0001 and .9999")}
 if (p2 > .9999 || p2 < .0001) {stop("p2 must be between .0001 and .9999")}
 za <- qnorm(1 - alpha)
 zb <- qnorm(1 - (1 - pow)/2)
 cov <- phi*sqrt(p1*p2*(1 - p1)*(1 - p2))
 es <- p1 - p2
 n <- ceiling(((p1*(1 - p1) + p2*(1 - p2) - 2*cov))*(za + zb)^2/(h - abs(es))^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 rownames(out) <- ""
 return(out)
}


#  size.supinf.prop.ps =======================================================
#' Sample size for a paired-samples superiority or inferiority test of 
#' proportions
#'
#'
#' @description
#' Computes the sample size required to perform a superiority or inferiority 
#' test for the difference in population proportions with desired power in a
#' paired-samples design. For a superiority test, specify the upper limit (h)
#' for the range of practical equivalence and specify values of p1 and p2 such
#' that p1 - p2 > h. For an inferiority test, specify the lower limit (-h) for 
#' the range of practical equivalence and specify values of p1 and p2 such  
#' that p1 - p2 > -h. This function sets the effect size equal to the 
#' difference in proportion planning values. Set the phi correlation planning 
#' value to the smallest value within a plausible range for a conservatively 
#' large sample size.
#'
#'
#' @param  alpha  alpha level for hypothesis test 
#' @param  pow    desired power
#' @param  p1     planning value of proportion for measurement 1
#' @param  p2     planning value of proportion for measurement 2
#' @param  phi    planning value of phi correlation
#' @param  h      lower or upper limit for range of practical equivalence
#'
#'
#' @return
#' Returns the required sample size
#'
#'
#' @examples
#' size.supinf.prop.ps(.05, .9, .35, .20, .45, .05)
#'
#' # Should return:
#' # Sample size
#' #         227
#'
#'
#' @importFrom stats qnorm
#' @export
size.supinf.prop.ps <- function(alpha, pow, p1, p2, phi, h) {
 if (phi > .999 || phi < -.999) {stop("phi must be between -.999 and .999")}
 if (p1 > .9999 || p1 < .0001) {stop("p1 must be between .0001 and .9999")}
 if (p2 > .9999 || p2 < .0001) {stop("p2 must be between .0001 and .9999")}
 za <- qnorm(1 - alpha/2)
 zb <- qnorm(pow)
 cov <- phi*sqrt(p1*p2*(1 - p1)*(1 - p2))
 es <- p1 - p2
 n <- ceiling(((p1*(1 - p1) + p2*(1 - p2) - 2*cov))*(za + zb)^2/(es - h)^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 rownames(out) <- ""
 return(out)
}

          
# ===================== Power for Planned Sample Size  ========================
#  power.prop ================================================================
#' Approximates the power of a 1-group proportion test for a planned sample
#' size
#'
#'
#' @description
#' Computes the approximate power of a one-sample proportion test for a
#' planned sample size. Set the proportion planning value to .5 for a 
#' conservatively low power estimate. The value of the effect size need
#' not be based on the proportion planning value.
#'
#'
#' @param  alpha  alpha level for hypothesis test 
#' @param  n      planned sample size
#' @param  p      planning value of proportion 
#' @param  es     planning value of proportion minus null hypothesis value
#'
#'
#' @return
#' Returns the approximate power of the test
#'
#'
#' @examples
#' power.prop(.05, 40, .5, .2)
#'
#' # Should return:
#' #     Power
#' # 0.7156044
#'
#'
#' @importFrom stats qnorm
#' @importFrom stats pnorm
#' @export
power.prop <- function(alpha, n, p, es) {
 if (p > .9999 || p < .0001) {stop("proportion must be between .0001 and .9999")}
 za <- qnorm(1 - alpha/2)
 z1 <- abs(es)/sqrt(p*(1 - p)/n) - za
 z2 <- abs(es)/sqrt(p*(1 - p)/n) + za
 pow1 <- pnorm(z1)
 pow2 <- 1 - pnorm(z2)
 pow <- pow1 + pow2
 out <- matrix(pow, nrow = 1, ncol = 1)
 colnames(out) <- "Power"
 rownames(out) <- ""
 return(out)
}

          
#  power.prop2 ================================================================
#' Approximates the power of a 2-group proportion test for planned sample 
#' sizes
#'
#'
#' @description
#' Computes the approximate power for a test of equal population proportions
#' in a 2-group design for the planned sample sizes. This function requires 
#' planning values for both proportions. Set the proportion planning values 
#' to .5 for a conservatively low power estimate. The planning value for the 
#' proportion difference could be set to the difference of the two proportion
#' planning values or it could be set to a minimally interesting effect size.
#'
#'
#' @param  alpha  alpha level for hypothesis test 
#' @param  n1     planned sample size for group 1
#' @param  n2     planned sample size for group 2
#' @param  p1     planning value of proportion for group 1
#' @param  p2     planning value of proportion for group 2
#' @param  es     planning value of proportion difference
#'
#'
#' @return
#' Returns the approximate power of the test
#'
#'
#' @examples
#' power.prop2(.05, 60, 40, .5, .5, .2)
#'
#' # Should return:
#' #     Power
#' # 0.4998959
#'
#'
#' @importFrom stats qnorm
#' @importFrom stats pnorm
#' @export
power.prop2 <- function(alpha, n1, n2, p1, p2, es) {
 if (p1 > .9999 || p1 < .0001) {stop("p1 must be between .0001 and .9999")}
 if (p2 > .9999 || p2 < .0001) {stop("p2 must be between .0001 and .9999")}
 za <- qnorm(1 - alpha/2)
 z1 <- abs(es)/sqrt(p1*(1 - p1)/n1 + p2*(1 - p2)/n2) - za
 z2 <- abs(es)/sqrt(p1*(1 - p1)/n1 + p2*(1 - p2)/n2) + za
 pow1 <- pnorm(z1)
 pow2 <- 1 - pnorm(z2)
 pow <- pow1 + pow2
 out <- matrix(pow, nrow = 1, ncol = 1)
 colnames(out) <- "Power"
 rownames(out) <- ""
 return(out)
}


#  power.prop.ps ==============================================================
#' Approximates the power of a paired-samples test of equal porportions for a
#' planned sample size
#'
#'                     
#' @description
#' Computes the approximate power of a test for equal population proportions in
#' a paired-samples design (the McNemar test). This function requires planning 
#' values for both proportions and a phi coefficient that describes the 
#' correlation between the two dichotomous measurements. The proportion planning 
#' values can set to .5 for a conservatively low power estimate. The planning 
#' value for the proportion difference (effect size) could be set to the
#' difference of the two proportion planning values or it could be set to a
#' minimally interesting effect size. Set the phi correlation planning value 
#' to the smallest value within a plausible range for a conservatively low power
#' estimate. 
#'
#'
#' @param  alpha  alpha level for hypothesis test 
#' @param  n      planned sample size
#' @param  p1     planning value of proportion for measurement 1
#' @param  p2     planning value of proportion for measurement 2
#' @param  phi    planning value of phi correlation
#' @param  es     planning value of proportion difference
#'
#'
#' @return
#' Returns the approximate power of the test
#'
#'
#' @examples
#' power.prop.ps(.05, 45, .5, .5, .4, .2)
#'
#' # Should return:
#' #     Power
#' # 0.6877704
#'
#'
#' @importFrom stats qnorm
#' @importFrom stats pnorm
#' @export
power.prop.ps <- function(alpha, n, p1, p2, phi, es) {
 if (phi > .999 || phi < -.999) {stop("phi must be between -.999 and .999")}
 if (p1 > .9999 || p1 < .0001) {stop("p1 must be between .0001 and .9999")}
 if (p2 > .9999 || p2 < .0001) {stop("p2 must be between .0001 and .9999")}
 za <- qnorm(1 - alpha/2)
 v <- 2*phi*sqrt(p1*(1 - p1)*p2*(1 - p2))
 z1 <- abs(es)/sqrt((p1*(1 - p1) + p2*(1 - p2) - v)/n) - za
 z2 <- abs(es)/sqrt((p1*(1 - p1) + p2*(1 - p2) - v)/n) + za
 pow1 <- pnorm(z1)
 pow2 <- 1 - pnorm(z2)
 pow <- pow1 + pow2
 out <- matrix(pow, nrow = 1, ncol = 1)
 colnames(out) <- "Power"
 rownames(out) <- ""
 return(out)
}

          
# ======================== Miscellaneous =====================================
#  iqv =======================================================================
#' Indices of qualitative variation 
#'
#'
#' @description
#' Computes the Shannon, Berger, and Simpson indices of qualitative variation.
#'
#'  
#' @param   f   vector of multinomial frequency counts
#'
#' 
#' @return 
#' Returns estimates of the Shannon, Berger, and Simpson indices
#' 
#' 
#' @examples
#' f <- c(10, 46, 15, 3)
#' iqv(f)
#'
#' # Should return:
#' #   Simpson    Berger   Shannon
#' # 0.7367908 0.5045045       0.7
#'  
#' 
#' @export  
iqv <- function(f) {
 n <- sum(f)
 p <- f/n
 k <- length(f)
 a <- k/(k - 1)
 iqv1 <- a*(1 - sum(p^2))
 iqv2 <- a*( 1 - max(p))
 iqv3 <- (-1)*sum(p*log(p))/log(k)
 out <- t(c(iqv1, iqv2, iqv3))
 colnames(out) <- c("Simpson", "Berger", "Shannon")
 rownames(out) <- ""
 return(out)
}


fix_imports <- function() {
  something <- Rdpack::append_to_Rd_list()
  res <- mathjaxr::preview_rd()
}
