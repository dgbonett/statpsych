#  ================== New functions for version 1.2 of statpsych ===================
#
#  ci.cramer  ======================================================================
#' Confidence interval for Cramer's V
#'
#'
#' @description
#' Computes a confidence interval for a population Cramer's V coefficient
#' of nominal association for an r x s contingency table and its approximate
#' standard error. The confidence interval is based on a noncentral chi-square 
#' distribution, and an approximate standard error is extracted from the
#' confidence interval.
#'
#'
#' @param  alpha    alpha value for 1-alpha confidence
#' @param  chisqr   Pearson chi-square test statistic for independence
#' @param  r        number of rows in contingency table
#' @param  c        number of columns in contengency table
#' @param  n        sample size
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Cramer's V - estimate of Cramer's V 
#' * SE - approximate standard error 
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
#' #      Cramer's V     SE     LL     UL
#' # [1,]     0.3099 0.0674 0.1888 0.4529
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
  du <- n*k - df
  nc <- seq(0, du, by = .001)
  p <- pchisq(chisqr, df, nc)
  k1 <- which(min(abs(p - alpha2)) == abs(p - alpha2))[[1]]
  dL <- nc[k1]
  LL <- sqrt((dL + df)/(n*k))
  k2 <- which(min(abs(p - alpha1)) == abs(p - alpha1))[[1]]
  dU <- nc[k2]
  UL <- sqrt((dU + df)/(n*k))
  se <- (UL - LL)/(2*z)
  out <- round(t(c(v, se, LL, UL)), 4)
  colnames(out) <- c("Cramer's V", "SE", "LL", "UL")
  return(out)
}


#  sim.ci.mean1 ===============================================================
#' Simulates confidence interval coverage probability for single mean
#'
#'                               
#' @description
#' Performs a computer simulation (20,000 Monte Carlo samples) of the 
#' confidence interval performance for a single mean. Sample data can 
#' be generated from five different population distributions. All 
#' distributions are scaled to have standard deviations of 1.0.
#'
#' @param   alpha  alpha level for 1-alpha confidence
#' @param   n      sample size
#' @param   dist   type of distribution (1, 2, 3, 4,or 5)
#' * 1 = Gaussian (skewness = 0 and excess kurtosis = 0) 
#' * 2 = platykurtic (skewness = 0 and excess kurtosis = -1.2)
#' * 3 = leptokurtic (skewness = 0 and excess kurtsois = 6)
#' * 4 = moderate skew (skewness = 1 and excess kurtosis = 1.5)
#' * 5 = large skew (skewness = 2 and excess kurtosis = 6)
#'  
#' 
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Coverage - probability of confidence interval including population mean  
#' * Lower Error - probability of lower limit greater than population mean
#' * Upper Error - probability of upper limit less than population mean
#' * Ave CI Width - average confidence interval width
#'
#'
#' @examples
#' sim.ci.mean1(.05, 40, 4)
#'
#' # Should return (within sampling error):
#' #      Coverage Lower Error Upper Error Ave CI Width
#' # [1,]  0.94722     0.01738      0.0354    0.6333067
#'
#'
#' @importFrom stats qt
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom stats rt
#' @importFrom stats rgamma
#' @export
sim.ci.mean1 <- function(alpha, n, dist) {
  tcrit <- qt(1 - alpha/2, n - 1)
  rep <- 20000
  if (dist == 1) {
    y <- matrix(rnorm(rep*n), nrow = rep)
    popmean <- 0
  } else if (dist == 2) {
    y <- matrix(3.464*runif(rep*n), nrow = rep)
    popmean <- 1.732
  } else if (dist == 3) {
    y <- matrix(.7745*rt(rep*n, 5), nrow = rep)
    popmean <- 0
  } else if (dist == 4) {
    y <- matrix(.5*rgamma(rep*n, 4), nrow = rep)
    popmean <- 2
  } else {
    y <- matrix(rgamma(rep*n, 1), nrow = rep)
    popmean <- 1
  }
  m <- rowMeans(y)
  var <- apply(y, 1, var)
  se <- sqrt(var/n)
  ll <- m - tcrit*se
  ul <- m + tcrit*se
  w <- ul - ll
  c1 <- as.integer(ll > popmean)
  c2 <- as.integer(ul < popmean)
  e1 <- sum(c1)/rep
  e2 <- sum(c2)/rep
  width <- mean(w)
  cov <- 1 - (e1 + e2)
  out <- t(c(cov, e1, e2, width))
  colnames(out) <- c("Coverage", "Lower Error", "Upper Error", "Ave CI Width")
  return(out)
}


#  sim.ci.mean2 ===============================================================
#' Simulates confidence interval coverage probability for a two-group mean 
#' difference
#'
#'                               
#' @description
#' Performs a computer simulation (20,000 Monte Carlo samples) of separate 
#' variance and pooled variance confidence interval performance for a mean  
#' difference in a two-group design. Sample data within each group can be 
#' generated from five different population distributions. All distributions
#' are scaled to have a standard deviation of 1.0 in group 1. 
#'
#' @param   alpha     alpha level for 1-alpha confidence
#' @param   n1        sample size in group 1
#' @param   n2        sample size in group 2
#' @param   sd.ratio  ratio of population standard deviations
#' @param   dist1     type of distribution in group 1 (1, 2, 3, 4, or 5)
#' @param   dist2     type of distribution in group 2 (1, 2, 3, 4, or 5)
#' * 1 = Gaussian (skewness = 0 and excess kurtosis = 0) 
#' * 2 = platykurtic (skewness = 0 and excess kurtosis = -1.2)
#' * 3 = leptokurtic (skewness = 0 and excess kurtsois = 6)
#' * 4 = moderate skew (skewness = 1 and excess kurtosis = 1.5)
#' * 5 = large skew (skewness = 2 and excess kurtosis = 6)
#'  
#' 
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Coverage - probability of confidence interval including population mean  
#' * Lower Error - probability of lower limit greater than population mean
#' * Upper Error - probability of upper limit less than population mean
#' * Ave CI Width - average confidence interval width
#'
#'
#' @examples
#' sim.ci.mean2(.05, 30, 25, 1.5, 4, 5)
#'
#' # Should return (within sampling error):
#' #                             Coverage Lower Error Upper Error Ave CI Width
#' # Equal Variances Assumed:      0.93986     0.04022     0.01992     1.344437
#' # Equal Variances Not Assumed:  0.94762     0.03862     0.01376     1.401305
#'
#'
#' @importFrom stats qt
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom stats rt
#' @importFrom stats rgamma
#' @export
sim.ci.mean2 <- function(alpha, n1, n2, sd.ratio, dist1, dist2) {
  df1 <- n1 + n2 - 2
  tcrit1 <- qt(1 - alpha/2, df1)
  rep <- 20000
  if (dist1 == 1) {
    y1 <- matrix(rnorm(rep*n1), nrow = rep)
    popmean1 <- 0
  } else if (dist1 == 2) {
    y1 <- matrix(3.464*runif(rep*n1), nrow = rep)
    popmean1 <- 1.732
  } else if (dist1 == 3) {
    y1 <- matrix(.7745*rt(rep*n1, 5), nrow = rep)
    popmean1 <- 0
  } else if (dist1 == 4) {
    y1 <- matrix(.5*rgamma(rep*n1, 4), nrow = rep)
    popmean1 <- 2
  } else {
    y1 <- matrix(rgamma(rep*n1, 1), nrow = rep)
    popmean1 <- 1
  }
  if (dist2 == 1) {
    y2 <- matrix(sd.ratio*rnorm(rep*n2), nrow = rep)
    popmean2 <- 0
  } else if (dist2 == 2) {
    y2 <- matrix(sd.ratio*3.464*runif(rep*n2), nrow = rep)
    popmean2 <- sd.ratio*1.732
  } else if (dist2 == 3) {
    y2 <- matrix(sd.ratio*.7745*rt(rep*n2, 5), nrow = rep)
    popmean2 <- 0
  } else if (dist2 == 4) {
    y2 <- matrix(sd.ratio*.5*rgamma(rep*n2, 4), nrow = rep)
    popmean2 <- sd.ratio*2
  } else {
    y2 <- matrix(sd.ratio*rgamma(rep*n2, 1), nrow = rep)
    popmean2<- sd.ratio
  }
  popdiff <- popmean1 - popmean2
  m1 <- rowMeans(y1)
  m2 <- rowMeans(y2)
  v1 <- apply(y1, 1, var)
  v2 <- apply(y2, 1, var)
  vp <- ((n1 - 1)*v1 + (n2 - 1)*v2)/df1
  se1 <- sqrt(vp/n1 + vp/n2)
  se2 <- sqrt(v1/n1 + v2/n2)
  df2 <- (se2^4)/(v1^2/(n1^3 - n1^2) + v2^2/(n2^3 - n2^2))
  tcrit2 <- qt(1 - alpha/2, df2)
  LL1 <- m1 - m2 - tcrit1*se1
  UL1 <- m1 - m2 + tcrit1*se1
  w1 <- UL1 - LL1
  LL2 <- m1 - m2 - tcrit2*se2
  UL2 <- m1 - m2 + tcrit2*se2
  w2 <- UL2 - LL2
  c11 <- as.integer(LL1 > popdiff)
  c12 <- as.integer(UL1 < popdiff)
  c21 <- as.integer(LL2 > popdiff)
  c22 <- as.integer(UL2 < popdiff)
  e11 <- sum(c11)/rep
  e12 <- sum(c12)/rep
  e21 <- sum(c21)/rep
  e22 <- sum(c22)/rep
  width1 <- mean(w1)
  width2 <- mean(w2)
  cov1 <- 1 - (e11 + e12)
  cov2 <- 1 - (e21 + e22)
  out1 <- t(c(cov1, e11, e12, width1))
  out2 <- t(c(cov2, e21, e22, width2))
  out <- rbind(out1, out2)
  colnames(out) <- c("Coverage", "Lower Error", "Upper Error", "Ave CI Width")
  rownames(out) <- c("Equal Variances Assumed:", "Equal Variances Not Assumed:")
  return(out)
}


#  sim.ci.mean.ps ===============================================================
#' Simulates confidence interval coverage probability for a paired-samples mean 
#' difference
#'
#'                               
#' @description
#' Performs a computer simulation (20,000 Monte Carlo samples) of confidence
#' interval performance for a mean difference in a paired-samples design. Sample 
#' data within each level of the within-subjects factor group can be generated
#' from bivariate population distributions with five different marginal
#' distributions. All distributions are scaled to have standard deviations of 
#' 1.0 at level 1. Bivariate random data with specified marginal skewness and
#' kurtosis are generated using the unonr function in the mnonr package. 
#'
#' @param   alpha     alpha level for 1-alpha confidence
#' @param   n         sample size 
#' @param   sd.ratio  ratio of population standard deviations
#' @param   cor       population correlation of paired observations
#' @param   dist1     type of distribution at level 1 (1, 2, 3, 4, or 5)
#' @param   dist2     type of distribution at level 2 (1, 2, 3, 4, or 5)
#' * 1 = Gaussian (skewness = 0 and excess kurtosis = 0) 
#' * 2 = platykurtic (skewness = 0 and excess kurtosis = -1.2)
#' * 3 = leptokurtic (skewness = 0 and excess kurtsois = 6)
#' * 4 = moderate skew (skewness = 1 and excess kurtosis = 1.5)
#' * 5 = large skew (skewness = 2 and excess kurtosis = 6)
#'  
#' 
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Coverage - probability of confidence interval including population mean  
#' * Lower Error - probability of lower limit greater than population mean
#' * Upper Error - probability of upper limit less than population mean
#' * Ave CI Width - average confidence interval width
#'
#'
#' @examples
#' library(mnonr)
#' sim.ci.mean.ps(.05, 30, 1.5, .7, 4, 5)
#'
#' # Should return (within sampling error):
#' #      Coverage Lower Error Upper Error Ave CI Width
#' # [1,]  0.93815     0.05125      0.0106    0.7778518
#'
#'
#' @importFrom stats qt
#' @importFrom mnonr unonr
#' @export
sim.ci.mean.ps <- function(alpha, n, sd.ratio, cor, dist1, dist2) {
  tcrit <- qt(1 - alpha/2, n - 1)
  rep <- 20000
  if (dist1 == 1) {
    skw1 <- 0; kur1 <- 0
  } else if (dist1 == 2) {
    skw1 <- 0; kur1 <- -1.2
  } else if (dist1 == 3) {
    skw1 <- 0; kur1 <- 6
  } else if (dist1 == 4) {
    skw1 <- .75; kur1 <- .86
  } else {
    skw1 <- 1.41; kur1 <- 3
  }
  if (dist2 == 1) {
    skw2 <- 0; kur2 <- 0
  } else if (dist2 == 2) {
    skw2 <- 0; kur2 <- -1.2
  } else if (dist2 == 3) {
    skw2 <- 0; kur2 <- 6
  } else if (dist2 == 4) {
    skw2 <- 1; kur2 <- 1.5
  } else {
    skw2 <- 2; kur2 <- 6
  }
  V <- matrix(c(1, cor*sd.ratio, cor*sd.ratio, sd.ratio^2), 2, 2)
  w <- 0; k <- 0; e1 <-0; e2 <- 0
  repeat {
    k <- k + 1
    y <- unonr(n, c(0, 0), V, skewness = c(skw1, skw2), kurtosis = c(kur1, kur2))
    q <- c(1, -1)
    d <- y%*%q
    m <- mean(d)
    se <- sqrt(var(d)/n)
    ll <- m - tcrit*se
    ul <- m + tcrit*se
    w0 <- ul - ll
    c1 <- as.integer(ll > 0)
    c2 <- as.integer(ul < 0)
    e1 <- e1 + c1
    e2 <- e2 + c2
    w <- w + w0
    if (k == rep) {break}
  }
  width <- w/rep
  cov <- 1 - (e1 + e2)/rep
  out <- t(c(cov, e1/rep, e2/rep, width))
  colnames(out) <- c("Coverage", "Lower Error", "Upper Error", "Ave CI Width")
  return(out)
}


#  sim.ci.median1 =============================================================
#' Simulates confidence interval coverage probability for single median
#'
#'                                      
#' @description
#' Performs a computer simulation (20,000 Monte Carlo samples) of the 
#' confidence interval performance for a single mean. Sample data can 
#' be generated from five different population distributions. All 
#' distributions are scaled to have standard deviations of 1.0.
#'
#' @param   alpha  alpha level for 1-alpha confidence
#' @param   n      sample size
#' @param   dist   type of distribution (1, 2, 3, 4, or 5) 
#' * 1 = Gaussian (skewness = 0 and excess kurtosis = 0) 
#' * 2 = platykurtic (skewness = 0 and excess kurtosis = -1.2)
#' * 3 = leptokurtic (skewness = 0 and excess kurtsois = 6)
#' * 4 = moderate skew (skewness = 1 and excess kurtosis = 1.5)
#' * 5 = large skew (skewness = 2 and excess kurtosis = 6)
#'  
#' 
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Coverage - probability of confidence interval including population mean  
#' * Lower Error - probability of lower limit greater than population mean
#' * Upper Error - probability of upper limit less than population mean
#' * Ave CI Width - average confidence interval width
#'
#'
#' @examples
#' sim.ci.median1(.05, 20, 5)
#'
#' # Should return (within sampling error):
#' #      Coverage Lower Error Upper Error Ave CI Width
#' # [1,]   0.9589      0.0216      0.0195    0.9735528
#'
#'
#' @importFrom stats qnorm
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom stats rt
#' @importFrom stats rgamma
#' @export
sim.ci.median1 <- function(alpha, n, dist) {
  zcrit <- qnorm(1 - alpha/2)
  rep <- 20000
  o <- round((n - zcrit*sqrt(n))/2)
  if (o < 1) {o = 1}
  k <- 0; w <- 0; e1 <- 0; e2 <- 0
  repeat {
    k <- k + 1
    if (dist == 1) {
      y <- rnorm(n)
      popmedian <- 0
      y <- sort(y)
      ll <- y[o]
      ul <- y[n - o + 1]
      w0 <- ul - ll
      c1 <- as.integer(ll > popmedian)
      c2 <- as.integer(ul < popmedian)
      e1 <- e1 + c1
      e2 <- e2 + c2
      w <- w + w0
    } else if (dist == 2) {
      y <- 3.464*runif(n)
      popmedian <- 1.732
      y <- sort(y)
      ll <- y[o]
      ul <- y[n - o + 1]
      w0 <- ul - ll
      c1 <- as.integer(ll > popmedian)
      c2 <- as.integer(ul < popmedian)
      e1 <- e1 + c1
      e2 <- e2 + c2
      w <- w + w0
    } else if (dist == 3) {
      y <- .7745*rt(n, 5)
      popmedian <- 0
      y <- sort(y)
      ll <- y[o]
      ul <- y[n - o + 1]
      w0 <- ul - ll
      c1 <- as.integer(ll > popmedian)
      c2 <- as.integer(ul < popmedian)
      e1 <- e1 + c1
      e2 <- e2 + c2
      w <- w + w0
    } else if (dist == 4) {
      y <- .5*rgamma(n, 4)
      popmedian <- 1.837
      y <- sort(y)
      ll <- y[o]
      ul <- y[n - o + 1]
      w0 <- ul - ll
      c1 <- as.integer(ll > popmedian)
      c2 <- as.integer(ul < popmedian)
      e1 <- e1 + c1
      e2 <- e2 + c2
      w <- w + w0
    } else {
      y <- rgamma(n, 1)
      popmedian <- 0.690
      y <- sort(y)
      ll <- y[o]
      ul <- y[n - o + 1]
      w0 <- ul - ll
      c1 <- as.integer(ll > popmedian)
      c2 <- as.integer(ul < popmedian)
      e1 <- e1 + c1
      e2 <- e2 + c2
      w <- w + w0
    }
    if (k == rep) {break}
  }
  width <- w/rep
  cov <- 1 - (e1 + e2)/rep
  out <- t(c(cov, e1/rep, e2/rep, width))
  colnames(out) <- c("Coverage", "Lower Error", "Upper Error", "Ave CI Width")
  return(out)
}


#  sim.ci.cor ===============================================================
#' Simulates confidence interval coverage probability for a Pearson
#' correlation
#'
#'                               
#' @description
#' Performs a computer simulation (20,000 Monte Carlo samples) of confidence
#' interval performance for a Pearson correlation. A bias adjustment is used
#' to reduce the bias of the Fisher transformed Pearson correlation. Sample 
#' data can be generated from bivariate population distributions with five 
#' different marginal distributions. All distributions are scaled to have 
#' standard deviations of 1.0 at level 1. Bivariate random data with specified
#' marginal skewness and kurtosis are generated using the unonr function in
#' the mnonr package. 
#'
#' @param   alpha     alpha level for 1-alpha confidence
#' @param   n         sample size 
#' @param   cor       population Pearson correlation
#' @param   dist1     type of distribution at level 1 (1, 2, 3, 4, or 5)
#' @param   dist2     type of distribution at level 2 (1, 2, 3, 4, or 5)
#' * 1 = Gaussian (skewness = 0 and excess kurtosis = 0) 
#' * 2 = platykurtic (skewness = 0 and excess kurtosis = -1.2)
#' * 3 = leptokurtic (skewness = 0 and excess kurtsois = 6)
#' * 4 = moderate skew (skewness = 1 and excess kurtosis = 1.5)
#' * 5 = large skew (skewness = 2 and excess kurtosis = 6)
#'  
#' 
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Coverage - probability of confidence interval including population mean  
#' * Lower Error - probability of lower limit greater than population mean
#' * Upper Error - probability of upper limit less than population mean
#' * Ave CI Width - average confidence interval width
#'
#'
#' @examples
#' sim.ci.cor(.05, 30, .7, 4, 5)
#'
#' # Should return (within sampling error):
#' #      Coverage Lower Error Upper Error Ave CI Width
#' # [1,]  0.93815     0.05125      0.0106    0.7778518
#'
#'
#' @importFrom stats qt
#' @importFrom mnonr unonr
#' @export
sim.ci.cor <- function(alpha, n, cor, dist1, dist2) {
  zcrit <- qnorm(1 - alpha/2)
  rep <- 20000
  rep <- 500
  if (dist1 == 1) {
    skw1 <- 0; kur1 <- 0
  } else if (dist1 == 2) {
    skw1 <- 0; kur1 <- -1.2
  } else if (dist1 == 3) {
    skw1 <- 0; kur1 <- 6
  } else if (dist1 == 4) {
    skw1 <- .75; kur1 <- .86
  } else {
    skw1 <- 1.41; kur1 <- 3
  }
  if (dist2 == 1) {
    skw2 <- 0; kur2 <- 0
  } else if (dist2 == 2) {
    skw2 <- 0; kur2 <- -1.2
  } else if (dist2 == 3) {
    skw2 <- 0; kur2 <- 6
  } else if (dist2 == 4) {
    skw2 <- 1; kur2 <- 1.5
  } else {
    skw2 <- 2; kur2 <- 6
  }
  V <- matrix(c(1, cor, cor, 1), 2, 2)
  w <- 0; k <- 0; e1 <-0; e2 <- 0
  repeat {
    k <- k + 1
    y <- unonr(n, c(0, 0), V, skewness = c(skw1, skw2), kurtosis = c(kur1, kur2))
    R <- cor(y)
    r <- R[1,2]
    zr <- log((1 + r)/(1 - r))/2 - r/(2*(n - 1))
    se.z <- sqrt(1/((n - 3)))
    ll0 <- zr - zcrit*se.z
    ul0 <- zr + zcrit*se.z
    ll <- (exp(2*ll0) - 1)/(exp(2*ll0) + 1)
    ul <- (exp(2*ul0) - 1)/(exp(2*ul0) + 1)
    w0 <- ul - ll
    c1 <- as.integer(ll > cor)
    c2 <- as.integer(ul < cor)
    e1 <- e1 + c1
    e2 <- e2 + c2
    w <- w + w0
    if (k == rep) {break}
  }
  width <- w/rep
  cov <- 1 - (e1 + e2)/rep
  out <- t(c(cov, e1/rep, e2/rep, width))
  colnames(out) <- c("Coverage", "Lower Error", "Upper Error", "Ave CI Width")
  return(out)
}


#  sim.ci.spear ===============================================================
#' Simulates confidence interval coverage probability for a Spearman
#' correlation
#'
#' @description
#' Performs a computer simulation (20,000 Monte Carlo samples) of confidence
#' interval performance for a Spearman correlation. Sample data can be 
#' generated from bivariate population distributions with five different
#' marginal distributions. All distributions are scaled to have standard 
#' deviations of 1.0 at level 1. Bivariate random data with specified marginal
#' skewness and kurtosis are generated using the unonr function in the mnonr
#' package. 
#'
#' @param   alpha     alpha level for 1-alpha confidence
#' @param   n         sample size 
#' @param   cor       population Spearman correlation
#' @param   dist1     type of distribution at level 1 (1, 2, 3, 4, or 5)
#' @param   dist2     type of distribution at level 2 (1, 2, 3, 4, or 5)
#' * 1 = Gaussian (skewness = 0 and excess kurtosis = 0) 
#' * 2 = platykurtic (skewness = 0 and excess kurtosis = -1.2)
#' * 3 = leptokurtic (skewness = 0 and excess kurtsois = 6)
#' * 4 = moderate skew (skewness = 1 and excess kurtosis = 1.5)
#' * 5 = large skew (skewness = 2 and excess kurtosis = 6)
#'  
#' 
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Coverage - probability of confidence interval including population mean  
#' * Lower Error - probability of lower limit greater than population mean
#' * Upper Error - probability of upper limit less than population mean
#' * Ave CI Width - average confidence interval width
#'
#'
#' @examples
#' sim.ci.spear(.05, 30, .7, 4, 5)
#'
#' # Should return (within sampling error):
#' #      Coverage Lower Error Upper Error Ave CI Width
#' # [1,]  0.96235     0.01255      0.0251    0.4257299
#'
#'
#' @importFrom stats qt
#' @importFrom mnonr unonr
#' @export
sim.ci.spear <- function(alpha, n, cor, dist1, dist2) {
  zcrit <- qnorm(1 - alpha/2)
  rep <- 20000
  if (dist1 == 1) {
    skw1 <- 0; kur1 <- 0
  } else if (dist1 == 2) {
    skw1 <- 0; kur1 <- -1.2
  } else if (dist1 == 3) {
    skw1 <- 0; kur1 <- 6
  } else if (dist1 == 4) {
    skw1 <- .75; kur1 <- .86
  } else {
    skw1 <- 1.41; kur1 <- 3
  }
  if (dist2 == 1) {
    skw2 <- 0; kur2 <- 0
  } else if (dist2 == 2) {
    skw2 <- 0; kur2 <- -1.2
  } else if (dist2 == 3) {
    skw2 <- 0; kur2 <- 6
  } else if (dist2 == 4) {
    skw2 <- 1; kur2 <- 1.5
  } else {
    skw2 <- 2; kur2 <- 6
  }
  V <- matrix(c(1, cor, cor, 1), 2, 2)
  y <- unonr(100000, c(0, 0), V, skewness = c(skw1, skw2), kurtosis = c(kur1, kur2))
  popR <- cor(y, method = "spearman")
  popspear <- popR[1,2]
  w <- 0; k <- 0; e1 <-0; e2 <- 0
  repeat {
    k <- k + 1
    y <- unonr(n, c(0, 0), V, skewness = c(skw1, skw2), kurtosis = c(kur1, kur2))
    R <- cor(y, method = "spearman")
    r <- R[1,2]
    zr <- log((1 + r)/(1 - r))/2 
    se.z <- sqrt((1 + r^2/2)*(1 - r^2)^2/(n - 3))
    ll0 <- zr - zcrit*se.z/(1 - r^2)
    ul0 <- zr + zcrit*se.z/(1 - r^2)
    ll <- (exp(2*ll0) - 1)/(exp(2*ll0) + 1)
    ul <- (exp(2*ul0) - 1)/(exp(2*ul0) + 1)
    w0 <- ul - ll
    c1 <- as.integer(ll > popspear)
    c2 <- as.integer(ul < popspear)
    e1 <- e1 + c1
    e2 <- e2 + c2
    w <- w + w0
    if (k == rep) {break}
  }
  width <- w/rep
  cov <- 1 - (e1 + e2)/rep
  out <- t(c(cov, e1/rep, e2/rep, width))
  colnames(out) <- c("Coverage", "Lower Error", "Upper Error", "Ave CI Width")
  return(out)
}


# ci.2x2.mean.mixed ===========================================================
#' Computes tests and confidence intervals of effects in a 2x2 mixed design 
#' for means
#'
#'
#' @description
#' Computes confidence intervals and p-values for the AB interaction effect, 
#' main effect of A, main efect of B, simple main effects of A, and simple main
#' effects of B in a 2x2 mixed factorial design with a quantitative response
#' variable where Factor A is a within-subjects factor, and Factor B is a 
#' between-subjects factor. A Satterthwaite adjustment to the degrees of 
#' freedom is used and equality of population variances is not assumed.
#'
#'
#' @param   alpha   alpha level for 1-alpha confidence
#' @param   y11     vector of scores at level 1 in group 1
#' @param   y12     vector of scores at level 2 in group 1
#' @param   y21     vector of scores at level 1 in group 2
#' @param   y22     vector of scores at level 2 in group 2
#'
#'
#' @return
#' Returns a 7-row matrix (one row per effect). The columns are:
#' * Estimate - estimate of effect
#' * SE - standard error 
#' * t - t test statistic 
#' * df - degrees of freedom
#' * p - p-value 
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @examples
#' y11 = c(18, 19, 20, 17, 20, 16)
#' y12 = c(19, 18, 19, 20, 17, 16)
#' y21 = c(19, 16, 16, 14, 16, 18)
#' y22 = c(16, 10, 12,  9, 13, 15)
#' ci.2x2.mean.mixed(.05, y11, y12, y21, y22)
#'
#' # Should return:
#' #            Estimate        SE         t       df            p         LL        UL
#' # AB:      -3.8333333 0.9803627 -3.910117 8.346534 0.0041247610 -6.0778198 -1.588847
#' # A:        2.0833333 0.4901814  4.250128 8.346534 0.0025414549  0.9610901  3.205577
#' # B:        3.7500000 1.0226599  3.666908 7.601289 0.0069250119  1.3700362  6.129964
#' # A at b1:  0.1666667 0.8333333  0.200000 5.000000 0.8493605140 -1.9754849  2.308818
#' # A at b2:  4.0000000 0.5163978  7.745967 5.000000 0.0005732451  2.6725572  5.327443
#' # B at a1:  1.8333333 0.9803627  1.870056 9.943850 0.0911668588 -0.3527241  4.019391
#' # B at a2:  5.6666667 1.2692955  4.464419 7.666363 0.0023323966  2.7173445  8.615989
#'
#'
#' @importFrom stats qt
#' @importFrom stats pt
#' @export
ci.2x2.mean.mixed <- function(alpha, y11, y12, y21, y22) {
  n1 <- length(y11)
  n2 <- length(y21)
  diff1 <- y11 - y12
  diff2 <- y21 - y22
  ave1 <- (y11 + y12)/2
  ave2 <- (y21 + y22)/2
  vd1 <- var(diff1)
  vd2 <- var(diff2)
  va1 <- var(ave1)
  va2 <- var(ave2)
  est1 <- mean(diff1) - mean(diff2)
  se1 <- sqrt(vd1/n1 + vd2/n2)
  df1 <- (se1^4)/(vd1^2/(n1^3 - n1^2) + vd2^2/(n2^3 - n2^2))
  tcrit1 <- qt(1 - alpha/2, df1)
  t1 <- est1/se1
  p1 <- 2*(1 - pt(abs(t1), df1))
  LL1 <- est1 - tcrit1*se1
  UL1 <- est1 + tcrit1*se1
  row1 <- c(est1, se1, t1, df1, p1, LL1, UL1)
  est2 <- (mean(diff1) + mean(diff2))/2
  se2 <- sqrt(vd1/n1 + vd2/n2)/2
  df2 <- (se2^4)/(vd1^2/((n1^3 - n1^2)*16) + vd2^2/((n2^3 - n2^2)*16))
  tcrit2 <- qt(1 - alpha/2, df2)
  t2 <- est2/se2
  p2 <- 2*(1 - pt(abs(t2), df2))
  LL2 <- est2 - tcrit2*se2
  UL2 <- est2 + tcrit2*se2
  row2 <- c(est2, se2, t2, df2, p2, LL2, UL2)
  est3 <- mean(ave1) - mean(ave2)
  se3 <- sqrt(va1/n1 + va2/n2)
  df3 <- (se3^4)/(va1^2/(n1^3 - n1^2) + va2^2/(n2^3 - n2^2))
  tcrit3 <- qt(1 - alpha/2, df3)
  t3 <- est3/se3
  p3 <- 2*(1 - pt(abs(t3), df3))
  LL3 <- est3 - tcrit3*se3
  UL3 <- est3 + tcrit3*se3
  row3 <- c(est3, se3, t3, df3, p3, LL3, UL3)
  est4 <- mean(diff1)
  se4 <- sqrt(vd1/n1)
  df4 <- n1 - 1
  tcrit4 <- qt(1 - alpha/2, df4)
  t4 <- est4/se4
  p4 <- 2*(1 - pt(abs(t4), df4))
  LL4 <- est4 - tcrit4*se4
  UL4 <- est4 + tcrit4*se4
  row4 <- c(est4, se4, t4, df4, p4, LL4, UL4)
  est5 <- mean(diff2)
  se5 <- sqrt(vd2/n2)
  df5 <- n2 - 1
  tcrit5 <- qt(1 - alpha/2, df5)
  t5 <- est5/se5
  p5 <- 2*(1 - pt(abs(t5), df5))
  LL5 <- est5 - tcrit5*se5
  UL5 <- est5 + tcrit5*se5
  row5 <- c(est5, se5, t5, df5, p5, LL5, UL5)
  est6 <- mean(y11) - mean(y21)
  se6 <- sqrt(var(y11)/n1 + var(y21)/n2)
  df6 <- (se6^4)/(var(y11)^2/(n1^3 - n1^2) + var(y21)^2/(n2^3 - n2^2))
  tcrit6 <- qt(1 - alpha/2, df6)
  t6 <- est6/se6
  p6 <- 2*(1 - pt(abs(t6), df6))
  LL6 <- est6 - tcrit6*se6
  UL6 <- est6 + tcrit6*se6
  row6 <- c(est6, se6, t6, df6, p6, LL6, UL6)
  est7 <- mean(y12) - mean(y22)
  se7 <- sqrt(var(y12)/n1 + var(y22)/n2)
  df7 <- (se7^4)/(var(y12)^2/(n1^3 - n1^2) + var(y22)^2/(n2^3 - n2^2))
  tcrit7 <- qt(1 - alpha/2, df7)
  t7 <- est7/se7
  p7 <- 2*(1 - pt(abs(t7), df7))
  LL7 <- est7 - tcrit7*se7
  UL7 <- est7 + tcrit7*se7
  row7 <- c(est7, se7, t7, df7, p7, LL7, UL7)
  out <- rbind(row1, row2, row3, row4, row5, row6, row7)
  rownames(out) <- c("AB:", "A:", "B:", "A at b1:", "A at b2:", "B at a1:", "B at a2:")
  colnames(out) = c("Estimate", "SE", "t", "df", "p", "LL", "UL")
  return(out)
}


# ci.2x2.mean.ws =============================================================
#' Computes tests and confidence intervals of effects in a 2x2 within-subjects 
#' design for means
#'
#'
#' @description
#' Computes confidence intervals and p-values for the AB interaction effect, 
#' main effect of A, main efect of B, simple main effects of A, and simple main
#' effects of B in a 2x2 within-subjects design with a quantitative response
#' variable. 
#'
#'
#' @param   alpha   alpha level for 1-alpha confidence
#' @param   y11     vector of scores at level 1 of A and level 1 of B
#' @param   y12     vector of scores at level 1 of A and level 2 of B
#' @param   y21     vector of scores at level 2 of A and level 1 of B
#' @param   y22     vector of scores at level 2 of A and level 2 of B
#'
#'
#' @return
#' Returns a 7-row matrix (one row per effect). The columns are:
#' * Estimate - estimate of effect
#' * SE - standard error 
#' * t - t test statistic 
#' * df - degrees of freedom
#' * p - p-value 
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @examples
#' y11 = c(1,2,3,4,5,7,7)
#' y12 = c(1,0,2,4,3,8,7)
#' y21 = c(4,5,6,7,8,9,8)
#' y22 = c(5,6,8,7,8,9,9)
#' ci.2x2.mean.ws(.05, y11, y12, y21, y22)
#'
#' # Should return:
#' #             Estimate        SE          t df            p          LL          UL
#' # AB:       1.28571429 0.5654449  2.2738102  6 0.0633355395 -0.09787945  2.66930802
#' # A:       -3.21428571 0.4862042 -6.6109784  6 0.0005765210 -4.40398462 -2.02458681
#' # B:       -0.07142857 0.2296107 -0.3110855  6 0.7662600658 -0.63326579  0.49040865
#' # A at b1: -2.57142857 0.2973809 -8.6469203  6 0.0001318413 -3.29909331 -1.84376383
#' # A at b2: -3.85714286 0.7377111 -5.2285275  6 0.0019599725 -5.66225692 -2.05202879
#' # B at a1:  0.57142857 0.4285714  1.3333333  6 0.2308094088 -0.47724794  1.62010508
#' # B at a2: -0.71428571 0.2857143 -2.5000000  6 0.0465282323 -1.41340339 -0.01516804
#'
#'
#' @importFrom stats qt
#' @importFrom stats pt
#' @export
ci.2x2.mean.ws <- function(alpha, y11, y12, y21, y22) {
  n <- length(y11)
  df <- n - 1
  t <- qt(1 - alpha/2, df)
  q1 <- c(1, -1, -1, 1)
  q2 <- c(.5, .5, -.5, -.5)
  q3 <- c(.5, -.5, .5, -.5)
  q4 <- c(1, 0, -1, 0)
  q5 <- c(0, 1, 0, -1)
  q6 <- c(1, -1, 0, 0)
  q7 <- c(0, 0, 1, -1)
  y <- cbind(y11, y12, y21, y22)
  est1 <- mean(q1%*%t(y))
  se1 <- sqrt(var(matrix(q1%*%t(y)))/n)
  t1 <- est1/se1
  p1 <- 2*(1 - pt(abs(t1), df))
  LL1 <- est1 - t*se1
  UL1 <- est1 + t*se1
  row1 <- c(est1, se1, t1, df, p1, LL1, UL1)
  est2 <- mean(q2%*%t(y))
  se2 <- sqrt(var(matrix(q2%*%t(y)))/n)
  t2 <- est2/se2
  p2 <- 2*(1 - pt(abs(t2), df))
  LL2 <- est2 - t*se2
  UL2 <- est2 + t*se2
  row2 <- c(est2, se2, t2, df, p2, LL2, UL2)
  est3 <- mean(q3%*%t(y))
  se3 <- sqrt(var(matrix(q3%*%t(y)))/n)
  t3 <- est3/se3
  p3 <- 2*(1 - pt(abs(t3), df))
  LL3 <- est3 - t*se3
  UL3 <- est3 + t*se3
  row3 <- c(est3, se3, t3, df, p3, LL3, UL3)
  est4 <- mean(q4%*%t(y))
  se4 <- sqrt(var(matrix(q4%*%t(y)))/n)
  t4 <- est4/se4
  p4 <- 2*(1 - pt(abs(t4), df))
  LL4 <- est4 - t*se4
  UL4 <- est4 + t*se4
  row4 <- c(est4, se4, t4, df, p4, LL4, UL4)
  est5 <- mean(q5%*%t(y))
  se5 <- sqrt(var(matrix(q5%*%t(y)))/n)
  t5 <- est5/se5
  p5 <- 2*(1 - pt(abs(t5), df))
  LL5 <- est5 - t*se5
  UL5 <- est5 + t*se5
  row5 <- c(est5, se5, t5, df, p5, LL5, UL5)
  est6 <- mean(q6%*%t(y))
  se6 <- sqrt(var(matrix(q6%*%t(y)))/n)
  t6 <- est6/se6
  p6 <- 2*(1 - pt(abs(t6), df))
  LL6 <- est6 - t*se6
  UL6 <- est6 + t*se6
  row6 <- c(est6, se6, t6, df, p6, LL6, UL6)
  est7 <- mean(q7%*%t(y))
  se7 <- sqrt(var(matrix(q7%*%t(y)))/n)
  t7 <- est7/se7
  p7 <- 2*(1 - pt(abs(t7), df))
  LL7 <- est7 - t*se7
  UL7 <- est7 + t*se7
  row7 <- c(est7, se7, t7, df, p7, LL7, UL7)
  out <- rbind(row1, row2, row3, row4, row5, row6, row7)
  rownames(out) <- c("AB:", "A:", "B:", "A at b1:", "A at b2:", "B at a1:", "B at a2:")
  colnames(out) = c("Estimate", "SE", "t", "df", "p", "LL", "UL")
  return(out)
}


# ci.2x2.mean.bs =============================================================
#' Computes tests and confidence intervals of effects in a 2x2 betwen-subjects 
#' design for means
#'
#'
#' @description
#' Computes confidence intervals and p-values for the AB interaction effect, 
#' main effect of A, main efect of B, simple main effects of A, and simple main
#' effects of B in a 2x2 between-subjects design with a quantitative response
#' variable. A Satterthwaite adjustment to the degrees of freedom is used and 
#' equality of population variances is not assumed.
#'
#'
#' @param   alpha   alpha level for 1-alpha confidence
#' @param   y11     vector of scores at level 1 of A and level 1 of B
#' @param   y12     vector of scores at level 1 of A and level 2 of B
#' @param   y21     vector of scores at level 2 of A and level 1 of B
#' @param   y22     vector of scores at level 2 of A and level 2 of B
#'
#'
#' @return
#' Returns a 7-row matrix (one row per effect). The columns are:
#' * Estimate - estimate of effect
#' * SE - standard error 
#' * t - t test statistic 
#' * df - degrees of freedom
#' * p - p-value 
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @examples
#' y11 = c(14, 15, 11, 7, 16, 12, 15, 16, 10, 9)
#' y12 = c(18, 24, 14, 18, 22, 21, 16, 17, 14, 13)
#' y21 = c(16, 11, 10, 17, 13, 18, 12, 16, 6, 15)
#' y22 = c(18, 17, 11, 9, 9, 13, 18, 15, 14, 11)
#' ci.2x2.mean.bs(.05, y11, y12, y21, y22)
#'
#' # Should return:
#' #          Estimate       SE           t       df           p         LL         UL
#' # AB:         -5.10 2.224860 -2.29227953 35.47894 0.027931810 -9.6145264 -0.5854736
#' # A:           1.65 1.112430  1.48323970 35.47894 0.146840430 -0.6072632  3.9072632
#' # B:          -2.65 1.112430 -2.38217285 35.47894 0.022698654 -4.9072632 -0.3927368
#' # A at b1:    -0.90 1.545244 -0.58243244 17.56296 0.567678242 -4.1522367  2.3522367
#' # A at b2:     4.20 1.600694  2.62386142 17.93761 0.017246053  0.8362274  7.5637726
#' # B at a1:    -5.20 1.536952 -3.38331916 17.61093 0.003393857 -8.4341379 -1.9658621
#' # B at a2:    -0.10 1.608657 -0.06216365 17.91650 0.951120753 -3.4807927  3.2807927
#'
#'
#' @importFrom stats qt
#' @importFrom stats pt
#' @export
ci.2x2.mean.bs <- function(alpha, y11, y12, y21, y22) {
  n11 <- length(y11)
  n12 <- length(y12)
  n21 <- length(y21)
  n22 <- length(y22)
  v1 <- c(1, -1, -1, 1)
  v2 <- c(.5, .5, -.5, -.5)
  v3 <- c(.5, -.5, .5, -.5)
  v4 <- c(1, 0, -1, 0)
  v5 <- c(0, 1, 0, -1)
  v6 <- c(1, -1, 0, 0)
  v7 <- c(0, 0, 1, -1)
  m11 <- mean(y11)
  m12 <- mean(y12)
  m21 <- mean(y21)
  m22 <- mean(y22)
  sd11 <- sd(y11)
  sd12 <- sd(y12)
  sd21 <- sd(y21)
  sd22 <- sd(y22)
  m <- c(m11, m12, m21, m22)
  sd <- c(sd11, sd12, sd21, sd22) 
  n <- c(n11, n12, n21, n22)
  var <- diag(sd^2)%*%(solve(diag(n)))
  est1 <- t(v1)%*%m 
  se1 <- sqrt(t(v1)%*%var%*%v1)
  t1 <- est1/se1
  df1 <- (se1^4)/sum(((v1^4)*(sd^4)/(n^2*(n - 1))))
  tcrit1 <- qt(1 - alpha/2, df1)
  p1 <- 2*(1 - pt(abs(t1), df1))
  LL1 <- est1 - tcrit1*se1
  UL1 <- est1 + tcrit1*se1
  row1 <- c(est1, se1, t1, df1, p1, LL1, UL1)
  est2 <- t(v2)%*%m 
  se2 <- sqrt(t(v2)%*%var%*%v2)
  t2 <- est2/se2
  df2 <- (se2^4)/sum(((v2^4)*(sd^4)/(n^2*(n - 1))))
  tcrit2 <- qt(1 - alpha/2, df2)
  p2 <- 2*(1 - pt(abs(t2), df2))
  LL2 <- est2 - tcrit2*se2
  UL2 <- est2 + tcrit2*se2
  row2 <- c(est2, se2, t2, df2, p2, LL2, UL2)
  est3 <- t(v3)%*%m 
  se3 <- sqrt(t(v3)%*%var%*%v3)
  t3 <- est3/se3
  df3 <- (se3^4)/sum(((v3^4)*(sd^4)/(n^2*(n - 1))))
  tcrit3 <- qt(1 - alpha/2, df3)
  p3 <- 2*(1 - pt(abs(t3), df3))
  LL3 <- est3 - tcrit3*se3
  UL3 <- est3 + tcrit3*se3
  row3 <- c(est3, se3, t3, df3, p3, LL3, UL3)
  est4 <- t(v4)%*%m 
  se4 <- sqrt(t(v4)%*%var%*%v4)
  t4 <- est4/se4
  df4 <- (se4^4)/sum(((v4^4)*(sd^4)/(n^2*(n - 1))))
  tcrit4 <- qt(1 - alpha/2, df4)
  p4 <- 2*(1 - pt(abs(t4), df4))
  LL4 <- est4 - tcrit4*se4
  UL4 <- est4 + tcrit4*se4
  row4 <- c(est4, se4, t4, df4, p4, LL4, UL4)
  est5 <- t(v5)%*%m 
  se5 <- sqrt(t(v5)%*%var%*%v5)
  t5 <- est5/se5
  df5 <- (se5^4)/sum(((v5^4)*(sd^4)/(n^2*(n - 1))))
  tcrit5 <- qt(1 - alpha/2, df5)
  p5 <- 2*(1 - pt(abs(t5), df5))
  LL5 <- est5 - tcrit5*se5
  UL5 <- est5 + tcrit5*se5
  row5 <- c(est5, se5, t5, df5, p5, LL5, UL5)
  est6 <- t(v6)%*%m 
  se6 <- sqrt(t(v6)%*%var%*%v6)
  t6 <- est6/se6
  df6 <- (se6^4)/sum(((v6^4)*(sd^4)/(n^2*(n - 1))))
  tcrit6 <- qt(1 - alpha/2, df6)
  p6 <- 2*(1 - pt(abs(t6), df6))
  LL6 <- est6 - tcrit6*se6
  UL6 <- est6 + tcrit6*se6
  row6 <- c(est6, se6, t6, df6, p6, LL6, UL6)
  est7 <- t(v7)%*%m 
  se7 <- sqrt(t(v7)%*%var%*%v7)
  t7 <- est7/se7
  df7 <- (se7^4)/sum(((v7^4)*(sd^4)/(n^2*(n - 1))))
  tcrit7 <- qt(1 - alpha/2, df7)
  p7 <- 2*(1 - pt(abs(t7), df7))
  LL7 <- est7 - tcrit7*se7
  UL7 <- est7 + tcrit7*se7
  row7 <- c(est7, se7, t7, df7, p7, LL7, UL7)
  out <- rbind(row1, row2, row3, row4, row5, row6, row7)
  rownames(out) <- c("AB:", "A:", "B:", "A at b1:", "A at b2:", "B at a1:", "B at a2:")
  colnames(out) = c("Estimate", "SE", "t", "df", "p", "LL", "UL")
  return(out)
}


# ci.2x2.prop.bs ===========================================================
#' Computes tests and confidence intervals of effects in a 2x2 between-
#' subjects design for proportions 
#'
#'
#' @description
#' Computes adjusted Wald confidence intervals and p-values for the AB 
#' interaction effect, main effect of A, main efect of B, simple main effects
#' of A, and simple main effects of B in a 2x2 between-subjects factorial 
#' design with a dichotomous response variable. The input vector of 
#' frequency counts is f11, f12, f21, f22, and the input vector of 
#' sample sizes is n11, n12, n21, n22 where the first subscript represents
#' the levels of Factor A and the second subscript represents the levels of
#' Factor B.
#'
#'
#' @param   alpha   alpha level for 1-alpha confidence
#' @param   f       vector of frequency counts of participants with attribute
#' @param   n       vector of sample sizes
#'
#'
#' @return
#' Returns a 7-row matrix (one row per effect). The columns are:
#' * Estimate - adjusted estimate of effect
#' * SE - standard error of estimate
#' * z - z test statistic for test of null hypothesis
#' * p - p-value 
#' * LL - lower limit of the adjusted Wald confidence interval
#' * UL - upper limit of the adjusted Wald confidence interval
#'
#'
#' @examples
#' f = c(15, 24, 28, 23)
#' n = c(50, 50, 50, 50)
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
#' Computes adjusted Wald confidence intervals and p-values for the AB 
#' interaction effect, main effect of A, main efect of B, simple main effects
#' of A, and simple main effects of B in a 2x2 mixed factorial design with a
#' dichotomous response variable where Factor A is a within-subjects factor 
#' and Factor B is a between-subjects factor. The 4x1 vector of frequency 
#' counts for Factor A within each group is  \[ f00, f01, f10, f11 \] where fij is 
#' the number of participants with a response of i = 0 or 1 at level 1 of 
#' Factor A and a response of j = 0 or 1 at level 2 of Factor A. 
#'
#'
#' @param   alpha   alpha level for 1-alpha confidence
#' @param   group1  2x2 contingency table for Factor A in group 1
#' @param   group2  2x2 contingency table for Factor A in group 2
#'
#'
#' @return
#' Returns a 7-row matrix (one row per effect). The columns are:
#' * Estimate - adjusted estimate of effect
#' * SE - standard error of estimate
#' * z - z test statistic 
#' * p - p-value
#' * LL - lower limit of the adjusted Wald confidence interval
#' * UL - upper limit of the adjusted Wald confidence interval
#'
#'
#' @examples
#' group1 = c(23, 42, 24, 11)
#' group2 = c(26, 27, 13, 34)
#' ci.2x2.prop.mixed (.05, group1, group2)
#'
#' # Should return:
#' #            Estimate         SE         z           p           LL        UL
#' # AB:      0.03960396 0.09991818 0.3963639 0.691836584 -0.156232072 0.2354400
#' # A:       0.15841584 0.04995909 3.1709113 0.001519615  0.060497825 0.2563339
#' # B:       0.09803922 0.04926649 1.9899778 0.046593381  0.001478675 0.1945998
#' # A at b1: 0.17647059 0.07893437 2.2356621 0.025373912  0.021762060 0.3311791
#' # A at b2: 0.13725490 0.06206620 2.2114274 0.027006257  0.015607377 0.2589024
#' # B at a1: 0.11764706 0.06842118 1.7194539 0.085531754 -0.016455982 0.2517501
#' # B at a2: 0.07843137 0.06913363 1.1344894 0.256589309 -0.057068054 0.2139308
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
  est1 <- (p1 - p2) - (p3 - p4)
  v1 <- (p1 + p2 - (p1 - p2)^2)/(n1 + 2)
  v2 <- (p3 + p4 - (p3 - p4)^2)/(n2 + 2)
  se1 <- sqrt(v1 + v2)
  z1 <- est1/se1
  pval1 <- 2*(1 - pnorm(abs(z1)))
  LL1 <- est1 - zcrit*se1
  UL1 <- est1 + zcrit*se1
  row1 <- c(est1, se1, z1, pval1, LL1, UL1)
  est2 <- ((p1 - p2) + (p3 - p4))/2
  se2 <- se1/2
  z2 <- est2/se2
  pval2 <- 2*(1 - pnorm(abs(z2)))
  LL2 <- est2 - zcrit*se2
  UL2 <- est2 + zcrit*se2
  row2 <- c(est2, se2, z2, pval2, LL2, UL2)
  p1 <- (2*f11 + f12 + f13 + 1)/(2*(n1 + 2))
  p2 <- (2*f21 + f22 + f23 + 1)/(2*(n2 + 2))
  est3 <- p1 - p2
  se3 <- sqrt(p1*(1 - p1)/(2*(n1 + 2)) + p2*(1 - p2)/(2*(n2 + 2)))
  z3 <- est3/se3
  pval3 <- 2*(1 - pnorm(abs(z3)))
  LL3 <- est3 - zcrit*se3
  UL3 <- est3 + zcrit*se3
  row3 <- c(est3, se3, z3, pval3, LL3, UL3)
  p1 <- (f12 + 1)/(n1 + 2)
  p2 <- (f13 + 1)/(n1 + 2)
  est4 <- p1 - p2
  se4 <- sqrt((p1 + p2 - (p1 - p2)^2)/(n1 + 2))
  z4 <- est4/se4
  pval4 <- 2*(1 - pnorm(abs(z4)))
  LL4 <- est4 - zcrit*se4
  UL4 <- est4 + zcrit*se4
  row4 <- c(est4, se4, z4, pval4, LL4, UL4)
  p1 <- (f22 + 1)/(n2 + 2)
  p2 <- (f23 + 1)/(n2 + 2)
  est5 <- p1 - p2
  se5 <- sqrt((p1 + p2 - (p1 - p2)^2)/(n2 + 2))
  z5 <- est5/se5
  pval5 <- 2*(1 - pnorm(abs(z5)))
  LL5 <- est5 - zcrit*se5
  UL5 <- est5 + zcrit*se5
  row5 <- c(est5, se5, z5, pval5, LL5, UL5)
  p1 <- (f11 + f12 + 1)/(n1 + 2)
  p2 <- (f21 + f22 + 1)/(n2 + 2)
  est6 <- p1 - p2
  se6 <- sqrt(p1*(1 - p1)/(n1 + 2) + p2*(1 - p2)/(n2 + 2))
  z6 <- est6/se6
  pval6 <- 2*(1 - pnorm(abs(z6)))
  LL6 <- est6 - zcrit*se6
  UL6 <- est6 + zcrit*se6
  row6 <- c(est6, se6, z6, pval6, LL6, UL6)
  p1 <- (f11 + f13 + 1)/(n1 + 2)
  p2 <- (f21 + f23 + 1)/(n2 + 2)
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
