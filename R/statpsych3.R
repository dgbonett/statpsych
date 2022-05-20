# ======================== File 3:  Confidence Intervals =====================
#  ci.prop1 ================================================================== 
#' Confidence interval for a single proportion
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


#  ci.pairs.prop1 ============================================================
#' Confidence intervals for pairwise proportion differences of a
#' polychotomous variable
#'
#'
#' @description
#' Computes adjusted Wald confidence intervals for pairwise proportion
#' differences of a polychotomous variable. These adjusted Wald confidence
#' intervals use the same method that is used to compare proportions in a 
#' paired-samples design.
#'
#'
#' @param   alpha   alpha level for 1-alpha confidence
#' @param   f       vector of multinomial frequency counts 
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - adjusted difference of proportions
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
#' ci.pairs.prop1(.05, f)
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
ci.pairs.prop1 <-function(alpha, f) {
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
#' @param   n1      sample size of group 1
#' @param   n2      sample size of group 2
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
#' #        Estimate         SE          LL        UL
#' # [1,] 0.09210526 0.04476077 0.004375769 0.1798348
#'
#'
#' @importFrom stats qnorm
#' @export
ci.prop2 <- function(alpha, f1, f2, n1, n2) {
 z <- qnorm(1 - alpha/2)
 p1 <- (f1 + 1)/(n1 + 2)
 p2 <- (f2 + 1)/(n2 + 2)
 se <- sqrt(p1*(1 - p1)/(n1 + 2) + p2*(1 - p2)/(n2 + 2))
 LL <- p1 - p2 - z*se
 UL <- p1 - p2 + z*se
 out <- t(c(p1-p2, se, LL, UL))
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
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
#' @param   n1      sample size of group 1
#' @param   n2      sample size of group 2
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - adjusted ratio of proportion ratio
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
#' #      Estimate       LL       UL
#' # [1,] 1.666667 1.017253 2.705025
#'
#'
#' @importFrom stats qnorm
#' @export
ci.ratio.prop2 <- function(alpha, f1, f2, n1, n2) {
 z <- qnorm(1 - alpha/2)
 p1 <- (f1 + 1/4)/(n1 + 7/4)
 p2 <- (f2 + 1/4)/(n2 + 7/4)
 v1 <- 1/(f1 + 1/4 + (f1 + 1/4)^2/(n1 - f1 + 3/2))
 v2 <- 1/(f2 + 1/4 + (f2 + 1/4)^2/(n2 - f2 + 3/2))
 se <- sqrt(v1 + v2)
 LL <- exp(log(p1/p2) - z*se)
 UL <- exp(log(p1/p2) + z*se)
 out <- t(c(p1/p2, LL, UL))
 colnames(out) <- c("Estimate", "LL", "UL")
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
#' @param   f       vector of frequency counts of participants with attribute
#' @param   n       vector of sample sizes
#' @param   v       vector of between-subjects contrast coefficients
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - adjusted estimate of proportion linear contrast
#' * SE - adjusted standard error
#' * z - z test statistic
#' * p - p-value
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
#' #       Estimate         SE        z           p         LL        UL
#' # [1,] 0.2119565 0.07602892 2.787841 0.005306059 0.06294259 0.3609705
#'
#'
#' @importFrom stats qnorm
#' @importFrom stats pnorm
#' @export
ci.lc.prop.bs <- function(alpha, f, n, v) {
 z <- qnorm(1 - alpha/2)
 m <- length(v) - length(which(v==0))
 p <- (f + 2/m)/(n + 4/m)
 est <- t(v)%*%p
 se <- sqrt(t(v)%*%diag(p*(1 - p))%*%solve(diag(n + 4/m))%*%v)
 zval <- est/se
 pval <- 2*(1 - pnorm(abs(zval)))
 LL <- est - z*se
 UL <- est + z*se
 CI <- c(LL, UL)
 out <- t(c(est, se, zval, pval, LL, UL))
 colnames(out) <- c("Estimate", "SE", "z", "p", "LL", "UL")
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
#' @param   f       vector of frequency counts of participants with attribute
#' @param   n       vector of sample sizes
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - adjusted estimate of proportion difference
#' * SE - adjusted standard error
#' * z - z test statistic
#' * p - p-value
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
#' Confidence interval for a slope of a quantitative between-subjects factor
#'
#'
#' @description
#' Computes an adjusted Wald confidence interval for the slope of a 
#' quantitative between-subjects factor using contrast coefficients. 
#'
#'
#' @param   alpha   alpha level for 1-alpha confidence
#' @param   f       vector of frequency counts of participants with attribute
#' @param   n       vector of sample sizes
#' @param   x       vector of quantitative factor values
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - adjusted slope estimate
#' * SE - adjusted standard error
#' * z - z test statistic
#' * p - p-value
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
#' #         Estimate          SE        z           p          LL         UL
#' # [1,] 0.007542293 0.002016793 3.739746 0.000184206 0.003589452 0.01149513
#'
#'
#' @importFrom stats qnorm
#' @importFrom stats pnorm
#' @export
ci.slope.prop.bs <- function(alpha, f, n, x) {
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
 LL <- slope - z*se
 UL <- slope + z*se
 out <- t(c(slope, se, t, pval, LL, UL))
 colnames(out) <- c("Estimate", "SE", "z", "p", "LL", "UL")
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
#' ci.prop.ps(.05, 12, 26, 4, 6)
#'
#' # Should return:
#' #       Estimate         SE        LL        UL
#' # [1,] 0.4583333 0.09448809 0.2548067 0.6251933
#'
#'
#' @importFrom stats qnorm
#' @export
ci.prop.ps <- function(alpha, f00, f01, f10, f11) {
 z <- qnorm(1 - alpha/2)
 n <- f00 + f01 + f10 + f11
 p1 <- (f01 + f11)/n
 p2 <- (f10 + f11)/n
 p01 <- (f01 + 1)/(n + 2)
 p10 <- (f10 + 1)/(n + 2)
 diff <- p1 - p2
 se <- sqrt(((p01 + p10) - (p01 - p10)^2)/(n + 2))
 LL <- p01 - p10 - z*se
 UL <- p01 - p10 + z*se
 out <- t(c(diff, se, LL, UL))
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
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
#' ci.ratio.prop.ps(.05, 12, 26, 4, 6)
#'
#' # Should return:
#' #      Estimate       LL       UL
#' # [1,]    2.375 1.537157 3.669518
#'
#'
#' @importFrom stats qnorm
#' @export
ci.ratio.prop.ps <- function(alpha, f00, f10, f01, f11) {
 z <- qnorm(1 - alpha/2)
 f1 <- f00 + f10
 f2 <- f00 + f01
 n0 <- f00 + f01 + f10
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
 LL <- exp(log(LL1) - log(UL2))
 UL <- exp(log(UL1) - log(LL2))
 out <- t(c(ratio, LL, UL))
 colnames(out) <- c("Estimate", "LL", "UL")
 return(out)
}


# ci.condslope.log ===========================================================
#' Confidence interval for conditional (simple) slopes in a logistic model
#'
#'
#' @description
#' Computes confidence intervals and test statistics for population 
#' conditional slopes (simple slopes) in a general linear model that
#' includes a predictor variable that is the product of a moderator 
#' variable and a predictor variable. Conditional slopes are computed 
#' at low and high values of the moderator variable. 
#'
#'
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  b1     sample coefficient for predictor variable
#' @param  b2     sample coefficient for product variable
#' @param  se1    standard error for predictor coefficient
#' @param  se2    standard error for product coefficient
#' @param  cov    sample covariance between predictor and product coefficients
#' @param  lo     low value of moderator variable 
#' @param  hi     high value of moderator variable 
#'
#'
#' @return 
#' Returns a 2-row matrix. The columns are:
#' * Estimate - estimated condition slope
#' * exp(Estimate) - estimated exponentiated condition slope
#' * z - z test statistic
#' * p - p-value
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
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
#' * Estimate - estimate of proportion ratio
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
#' #      Estimate       LL       UL
#' # [1,] 2.044451 1.133267 3.688254
#'
#'
#' @importFrom stats qnorm
#' @export
ci.oddsratio <- function(alpha, f00, f01, f10, f11) {
 z <- qnorm(1 - alpha/2)
 or <- (f11 + .5)*(f00 + .5)/((f01 + .5)*(f10 + .5))
 se.lor <- sqrt(1/(f00 + .5) + 1/(f01 + .5) + 1/(f10 + .5) + 1/(f11 + .5))
 LL <- exp(log(or) - z*se.lor)
 UL <- exp(log(or) + z*se.lor)
 out <- t(c(or, LL, UL))
 colnames(out) <- c("Estimate", "LL", "UL")
 return(out)
}


#  ci.yule ==================================================================== 
#' Confidence interval for Yule's Q
#'
#'
#' @description
#' Computes a confidence interval for Yule's Q measure of association using a
#' transformation of a confidence interval for an odds ratio with .5 added to
#' each cell frequency. This function requires the frequency counts from a 
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
#' * Estimate - estimate of Yule's Q
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @examples
#' ci.yule(.05, 229, 28, 96, 24)
#'
#' # Should return:
#' #      Estimate         LL       UL
#' # [1,] 0.343067 0.06247099 0.573402
#'
#'
#' @importFrom stats qnorm
#' @export
ci.yule <- function(alpha, f00, f01, f10, f11) {
 z <- qnorm(1 - alpha/2)
 or <- (f11 + .5)*(f00 + .5)/((f01 + .5)*(f10 + .5))
 se.lor <- sqrt(1/(f00 + .5) + 1/(f01 + .5) + 1/(f10 + .5) + 1/(f11 + .5))
 LLor <- exp(log(or) - z*se.lor)
 ULor <- exp(log(or) + z*se.lor)
 Q <- (or - 1)/(or + 1)
 LL <- (LLor - 1)/(LLor + 1)
 UL <- (ULor - 1)/(ULor + 1)
 out <- t(c(Q, LL, UL))
 colnames(out) <- c("Estimate", "LL", "UL")
 return(out)
}


#  ci.phi ====================================================================
#' Confidence interval for a phi correlation
#'
#'
#' @description
#' Computes a confidence interval for a phi correlation. This function requires 
#' the frequency counts from a 2 x 2 contingency table for two dichotomous 
#' variables. This measure of association assumes that both dichotomous  
#' variables are naturally dichotomous.
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
#' * Estimate - estimate of phi coefficient
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
#' #       Estimate         SE         LL        UL
#' # [1,] 0.1229976 0.05746271 0.01037273 0.2356224
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
 LL <- phi - z*se
 UL <- phi + z*se
 out <- t(c(phi, se, LL, UL))
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
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
#' @param   n1     sample size of group 1
#' @param   n2     sample size of group 2
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
#' #       Estimate         SE        LL       UL
#' # [1,] 0.4145733 0.07551281 0.2508866 0.546141
#'
#'
#' @importFrom stats qnorm
#' @export
ci.biphi <- function(alpha, f1, f2, n1, n2) {
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
 LL <- LL1/sqrt(LL1^2 + 2.89/(p1*p2))
 UL <- UL1/sqrt(UL1^2 + 2.89/(p1*p2))
 out <- t(c(biphi, se.biphi, LL, UL))
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
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
#' dichotomous. 
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
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @examples
#' ci.tetra(.05, 46, 15, 54, 85)
#'
#' # Should return:
#' #       Estimate        LL        UL
#' # [1,] 0.5135167 0.3102345 0.6748546
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
 LL <- cos(3.14159/(1 + LL1^c))
 UL <- cos(3.14159/(1 + UL1^c))
 out <- t(c(tetra, LL, UL))
 colnames(out) <- c("Estimate", "LL", "UL")
 return(out)
}


#  ci.kappa ================================================================== 
#' Confidence interval for a kappa reliability 
#'
#'
#' @description
#' Computes a confidence interval for the intraclass kappa coefficient and
#' Cohen's kappa coefficient for two dichotomous ratings. Both measures
#' are intraclass reliability coefficients.
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
#' Computes a confidence interval for a G-index of agreement between two
#' polychotomous ratings. This function requires the number of objects that
#' were given the same rating by both raters.
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
#' * Estimate - estimate of G-index of agreement 
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @examples
#' ci.agree(.05, 100, 80, 4)
#'
#' # Should return:
#' #       Estimate         SE        LL        UL
#' # [1,] 0.7333333 0.05333333 0.6132949 0.8226025
#'
#'
#' @importFrom stats qnorm
#' @export
ci.agree <- function(alpha, n, f, k) {
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
#' unobserved cell frequency.
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
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @examples
#' ci.popsize(.05, 794, 710, 741)
#'
#' # Should return:
#' #      Estimate   LL   UL
#' # [1,]     2908 2818 3012
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
 LL <- round(n0 + exp(log(f11) - z*se))
 UL <- round(n0 + exp(log(f11) + z*se))
 out <- t(c(N, LL, UL))
 colnames(out) <- c("Estimate", "LL", "UL")
 return(out)
}


# ======================== File 3: Hypothesis Tests ==========================
# test.prop1 =================================================================
#' Hypothesis test for a single proportion 
#'
#'
#' @description
#' Computes a continuity-corrected z test for a single proportion in a 
#' 1-group design.
#'
#'
#' @param   f    number of participants who have the attribute
#' @param   n    sample size
#' @param   h    population proportion under null hypothesis
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - ML estimate of proportion 
#' * z - z test statistic
#' * p - p-value
#'
#'
#' @references
#' \insertRef{Snedecor1980}{statpsych}
#'
#'
#' @examples
#' test.prop1(9, 20, .2)
#'
#' # Should return:
#' #      Estimate        z          p
#' # [1,]     0.45 2.515576 0.01188379
#'
#'
#' @importFrom stats pnorm
#' @export
test.prop1 <- function(f, n, h) {
 p <- f/n
 se <- sqrt(h*(1 - h)/n)
 z <- (abs(p - h) - 1/(2*n))/se
 pval <- 2*(1 - pnorm(abs(z)))
 out <- t(c(p, z, pval))
 colnames(out) <- c("Estimate", "z", "p")
 return(out)
}


# test.prop2 =================================================================
#' Hypothesis test for a 2-group proportion difference
#'
#'
#' @description
#' Computes a continuity-corrected z test for a difference of proportions in  
#' in a 2-group design.
#'
#'
#' @param  f1      number of group 1 participants who have the attribute
#' @param  f2      number of group 2 participants who have the attribute
#' @param  n1      sample size in group 1
#' @param  n2      sample size in group 2
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - ML estimate of proportion difference
#' * z - z test statistic
#' * p - p-value
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
#' #  Estimate        z           p
#' #      -0.3 2.899726 0.003734895
#'
#'
#' @importFrom stats pnorm
#' @export
test.prop2 <- function(f1, f2, n1, n2) {
 p1 <- f1/n1
 p2 <- f2/n2
 diff <- p1 - p2
 p0 <- (f1 + f2)/(n1 + n2)
 se <- sqrt(p0*(1 - p0)/n1 + p0*(1 - p0)/n2)
 z <- (abs(p1 - p2) - 1/(2*n1) - 1/(2*n2))/se
 pval <- 2*(1 - pnorm(abs(z)))
 out <- t(c(diff, z, pval))
 colnames(out) <- c("Estimate", "z", "p")
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
#' Chi-square - chi-square test statistic
#' df - degrees of freedom 
#' p - p-value
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
#' #      Chi-square df            p
#' # [1,]   17.41071  2 0.0001656958
#'
#'
#' @importFrom stats pchisq
#' @export
test.prop.bs <- function(f, n) {
 p0 <- sum(f)/sum(n)
 df <- length(f) - 1
 a <- 1/(p0*(1 - p0))
 p <- f/n
 r <- (p - p0)^2
 chi <- a*sum(n*r)
 pval <- 1 - pchisq(chi, df)
 out <- t(c(chi, df, pval))
 colnames(out) <- c("Chi-square", "df", "p")
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
#' Estimate - ML estimate of proportion difference
#' z - z test statistic
#' p - p-value
#'
#'
#' @examples
#' test.prop.ps(156, 96, 68, 80)
#'
#' # Should return:
#' #      Estimate        z          p
#' # [1,]     0.07 2.108346 0.03500109
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
 return(out)
}


# ================= File 3: Sample Size for Desired Precision ================
#  size.ci.prop1 
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
#' size.ci.prop1(.05, .4, .2)
#'
#' # Should return:
#' #      Sample size
#' # [1,]          93
#'
#'
#' @importFrom stats qnorm
#' @export                 
size.ci.prop1 <- function(alpha, p, w) {
 z <- qnorm(1 - alpha/2)
 n <- ceiling(4*p*(1 - p)*(z/w)^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
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
#' Returns the required sample size per group
#'
#'
#' @examples
#' size.ci.prop2(.05, .4, .2, .15)
#'
#' # Should return:
#' #      Sample size per group
#' # [1,]                   274
#'
#'
#' @importFrom stats qnorm
#' @export                 
size.ci.prop2 <- function(alpha, p1, p2, w) {
 z <- qnorm(1 - alpha/2)
 n <- ceiling(4*(p1*(1 - p1) + p2*(1 - p2))*(z/w)^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size per group"
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
#' Returns the required sample size per group
#'
#'
#' @examples
#' size.ci.ratio.prop2(.05, .2, .1, 2)
#'
#' # Should return:
#' #      Sample size per group
#' # [1,]                   416
#'
#'
#' @importFrom stats qnorm
#' @export   
size.ci.ratio.prop2 <- function(alpha, p1, p2, r) {
 z <- qnorm(1 - alpha/2)
 n <- ceiling(4*((1 - p1)/p1 + (1 - p2)/p2)*(z/log(r))^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size per group"
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
#' Returns the required sample size per group
#'
#'
#' @examples
#' p <- c(.25, .30, .50, .50)
#' v <- c(.5, .5, -.5, -.5)
#' size.ci.lc.prop.bs(.05, p, .2, v)
#'
#' # Should return:
#' #      Sample size per group
#' # [1,]                    87
#'
#'
#' @importFrom stats qnorm
#' @export
size.ci.lc.prop.bs <- function(alpha, p, w, v) {
 z <- qnorm(1 - alpha/2)
 n <- ceiling((4*t(v)%*%diag(p*(1 - p))%*%v)*(z/w)^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size per group"
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
#' @param  p1     planning value of proportion for group 1
#' @param  p2     planning value of proportion for group 2
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
#' #      Sample size
#' # [1,]         118
#'
#'
#' @importFrom stats qnorm
#' @export  
size.ci.prop.ps <- function(alpha, p1, p2, phi, w) {
 z <- qnorm(1 - alpha/2)
 cov <- phi*sqrt(p1*p2*(1 - p1)*(1 - p2))
 n <- ceiling((4*(p1*(1 - p1) + p2*(1 - p2) - 2*cov))*(z/w)^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 return(out)
}


#  size.ci.ratio.prop.ps ======================================================
#' Sample size for a paired-samples proportion ratio confidence interval  
#'
#'
#' @description
#' Computes the sample size required to estimate a ratio of proportions 
#' with desired confidence interval precision in a paired-samples design.
#' Set the phi planning value to the smallest value within a plausible range
#' for a conservatively large sample size. 
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
#' #      Sample size
#' # [1,]          67
#'
#'
#' @importFrom stats qnorm
#' @export  
size.ci.ratio.prop.ps <- function(alpha, p1, p2, phi, r) {
 z <- qnorm(1 - alpha/2)
 cov <- phi*sqrt((1 - p1)*(1 - p2)/(p1*p2))
 n <- ceiling(4*((1 - p1)/p1 + (1 - p2)/p2 - 2*cov)*(z/log(r))^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
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
#' #      Sample size
#' # [1,]         139
#'
#'
#' @importFrom stats qnorm
#' @export
size.ci.agree <- function(alpha, G, w) {
 z <- qnorm(1 - alpha/2)
 n <- ceiling(4*(1 - G^2)*(z/w)^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 return(out)
}


# ===================== File 3: Sample Size for Desired Power ================
#  size.test.prop1 =========================================================== 
#' Sample size for a test of a single proportion 
#'
#' @description
#' Computes the sample size required to test a single population proportion 
#' with desired power in a 1-group design. Set the proportion planning value 
#' to .5 for a conservatively large sample size.
#'
#'
#' @param  alpha  alpha level for hypothesis test 
#' @param  pow    desired power
#' @param  p      planning value of proportion 
#' @param  es     planning value of proportion minus null hypothesis value
#'
#'
#' @return
#' Returns the required sample size
#'
#'
#' @examples
#' size.test.prop1(.05, .9, .5, .2)
#'
#' # Should return:
#' #      Sample size
#' # [1,]          66
#'
#'
#' @importFrom stats qnorm
#' @export
size.test.prop1 <- function(alpha, pow, p, es) {
 za <- qnorm(1 - alpha/2)
 zb <- qnorm(pow)
 n <- ceiling(p*(1 - p)*(za + zb)^2/es^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 return(out)
}


#  size.test.prop2 ===========================================================
#' Sample size for a test of a 2-group proportion difference
#'
#'
#' @description
#' Computes the sample size in each group required to test a difference in 
#' population proportions with desired power in a 2-group design. This 
#' function requires planning values for both proportions. Set the proportion 
#' planning values to .5 for a conservatively large sample size. The
#' planning value for the proportion difference could be set equal to 
#' the difference of the two proportion planning values or it could be set
#' equal to a minimally interesting effect size.
#'
#'
#' @param  alpha  alpha level for hypothesis test 
#' @param  pow    desired power
#' @param  p1     planning value of proportion for group 1
#' @param  p2     planning value of proportion for group 2
#' @param  es     planning value of proportion difference
#'
#'
#' @return
#' Returns the required sample size per group
#'
#'
#' @examples
#' size.test.prop2(.05, .8, .2, .4, .2)
#'
#' # Should return:
#' #      Sample size per group
#' # [1,]                    79
#'
#'
#' @importFrom stats qnorm
#' @export
size.test.prop2 <- function(alpha, pow, p1, p2, es) {
 za <- qnorm(1 - alpha/2)
 zb <- qnorm(pow)
 n <- ceiling((p1*(1 - p1) + p2*(1 - p2))*(za + zb)^2/es^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size per group"
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
#' large sample size. The planning value for the linear contrast of proportions 
#' could be set equal to the linear contrast of proportion planning values or it
#' could be set equal to a minimally interesting effect size.
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
#' Returns the required sample size per group
#'
#'
#' @examples
#' p <- c(.25, .30, .50, .50)
#' v <- c(.5, .5, -.5, -.5)
#' size.test.lc.prop.bs(.05, .9, p, .15, v)
#'
#' # Should return:
#' #      Sample size per group
#' # [1,]                   105
#'
#'
#' @importFrom stats qnorm
#' @export
size.test.lc.prop.bs <- function(alpha, pow, p, es, v) {
 za <- qnorm(1 - alpha/2)
 zb <- qnorm(pow)
 n <- ceiling((t(v)%*%diag(p*(1 - p))%*%v)*(za + zb)^2/es^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size pe group"
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
#' Returns the required sample size per group
#'
#'
#' @examples
#' size.equiv.prop2(.1, .8, .30, .35, .15)
#'
#' # Should return:
#' #      Sample size per group
#' # [1,]                   288
#'
#'
#' @importFrom stats qnorm
#' @export
size.equiv.prop2 <- function(alpha, pow, p1, p2, h) {
 za <- qnorm(1 - alpha)
 zb <- qnorm(1 - (1 - pow)/2)
 es <- p1 - p2
 n <- ceiling((p1*(1 - p1) + p2*(1 - p2))*(za + zb)^2/(h - abs(es))^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
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
#' Returns the required sample size per group
#'
#'
#' @examples
#' size.supinf.prop2(.05, .9, .35, .20, .05)
#'
#' # Should return:
#' #      Sample size per group
#' # [1,]                   408
#'
#'
#' @importFrom stats qnorm
#' @export
size.supinf.prop2 <- function(alpha, pow, p1, p2, h) {
 za <- qnorm(1 - alpha/2)
 zb <- qnorm(pow)
 es <- p1 - p2
 n <- ceiling((p1*(1 - p1) + p2*(1 - p2))*(za + zb)^2/(es - h)^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
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
#' proportion planning values can set to .5 for a conservatively large sample 
#' size. The planning value for the proportion difference could be set equal 
#' to the difference of the two proportion planning values or it could be set
#' equal to a minimally interesting effect size. Set the phi correlation 
#' planning value to the smallest value within a plausible range for a 
#' conservatively large sample size.
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
#' #      Sample size
#' # [1,]         177
#'
#'
#' @importFrom stats qnorm
#' @export
size.test.prop.ps <- function(alpha, pow, p1, p2, phi, es) {
 za <- qnorm(1 - alpha/2)
 zb <- qnorm(pow)
 cov <- phi*sqrt(p1*p2*(1 - p1)*(1 - p2))
 n <- ceiling(((p1*(1 - p1) + p2*(1 - p2) - 2*cov))*((za + zb)/es)^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 return(out)
}


#  size.equiv.prop.ps ========================================================
#' Sample size for a paired-samples proportion equivalence test
#'
#'
#' @description
#' Computes the sample size in each group (assuming equal sample sizes) required
#' to perform an equivalence test for the difference in population proportions 
#' with desired power in a paired-samples design. The value of h specifies a 
#' range of practical equivalence, -h to h, for the difference in population 
#' proportions. The absolute difference in the proportion planning values must
#' be less than h.  Equivalence tests often require a very large sample size. 
#' Equivalence tests usually use 2 x alpha rather than alpha (e.g., use 
#' alpha = .10 rather alpha = .05). This function sets the effect size equal to
#' the difference in proportion planning values. Set the phi correlation planning
#' value to the smallest value within a plausible range for a conservatively 
#' large sample size.
#'
#'
#' @param  alpha  alpha level for hypothesis test 
#' @param  pow    desired power
#' @param  p1     planning value of proportion for group 1
#' @param  p2     planning value of proportion for group 2
#' @param  phi    planning value of phi coefficient
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
#' #      Sample size
#' # [1,]         173
#'
#'
#' @importFrom stats qnorm
#' @export
size.equiv.prop.ps <- function(alpha, pow, p1, p2, phi, h) {
 za <- qnorm(1 - alpha)
 zb <- qnorm(1 - (1 - pow)/2)
 cov <- phi*sqrt(p1*p2*(1 - p1)*(1 - p2))
 es <- p1 - p2
 n <- ceiling(((p1*(1 - p1) + p2*(1 - p2) - 2*cov))*(za + zb)^2/(h - abs(es))^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
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
#' #      Sample size
#' # [1,]         227
#'
#'
#' @importFrom stats qnorm
#' @export
size.supinf.prop.ps <- function(alpha, pow, p1, p2, phi, h) {
 za <- qnorm(1 - alpha/2)
 zb <- qnorm(pow)
 cov <- phi*sqrt(p1*p2*(1 - p1)*(1 - p2))
 es <- p1 - p2
 n <- ceiling(((p1*(1 - p1) + p2*(1 - p2) - 2*cov))*(za + zb)^2/(es - h)^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 return(out)
}


# ======================== File 3: Miscellaneous =============================
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
#' Returns estimates of the Shannon, Berger, and Simpson qualitative indices
#' 
#' 
#' @examples
#' f <- c(10, 46, 15, 3)
#' iqv(f)
#' # Should return:
#' #        Simpson    Berger   Shannon
#' # [1,] 0.7367908 0.5045045       0.7
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
 return(out)
}
