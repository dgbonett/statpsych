# ======================= File 1: Confidence Intervals =======================
#  ci.mean1  =================================================================
#' Confidence interval for a single mean
#'
#'
#' @description
#' Computes a confidence interval for a population mean using the sample 
#' mean, sample standard deviation, and sample size. Use the t.test function
#' for raw data input.
#'
#'  
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  m	  sample mean 
#' @param  sd	  sample standard deviation
#' @param  n	  sample size
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimated mean
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @references
#' \insertRef{Snedecor1980}{statpsych}
#'
#'
#' @examples
#' ci.mean1(.05, 24.5, 3.65, 40)
#'
#' # Should return:
#' #        Estimate        SE       LL       UL
#' # [1,]       24.5 0.5771157 23.33267 25.66733
#'  
#' 
#' @importFrom stats qt
#' @export
ci.mean1 <- function(alpha, m, sd, n) {
 df <- n - 1
 tcrit <- qt(1 - alpha/2, df)
 se <- sd/sqrt(n)
 ll <- m - tcrit*se
 ul <- m + tcrit*se
 out <- t(c(m, se, ll, ul))
 colnames(out) <- c("Estimate", "SE",  "LL", "UL")
 return(out)
}


#  ci.stdmean1  ==============================================================
#' Confidence interval for a single standardized mean
#'
#'
#' @description
#' Computes a confidence interval for a population standardized mean 
#' difference from a hypothesized value. If the hypothesized value is set
#' to 0, the reciprocals of the confidence interval endpoints gives a 
#' confidence interval for the coefficient of variation.
#'
#'  
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  m	  sample mean 
#' @param  sd	  sample standard deviation
#' @param  n	  sample size
#' @param  h      hypothesized value
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Estimate - bias adjusted standardized mean difference
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @references
#' \insertRef{Bonett2008}{statpsych}
#'
#'
#' @examples
#' ci.stdmean1(.05, 24.5, 3.65, 40, 20)
#'
#' # Should return:
#' #      Estimate        SE        LL       UL
#' # [1,] 1.209015 0.2124335 0.8165146 1.649239
#'  
#' 
#' @importFrom stats qt
#' @export
ci.stdmean1 <- function(alpha, m, sd, n, h) {
 z <- qnorm(1 - alpha/2)
 df <- n - 1
 adj <- 1 - 3/(4*df - 1)
 est <- (m - h)/sd
 estu <- adj*est
 se <- sqrt(est^2/(2*df) + 1/df)
 ll <- est - z*se
 ul <- est + z*se
 out <- t(c(estu, se, ll, ul))
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 return(out)
}
	

#  ci.mean2 ==================================================================
#' Confidence interval for a 2-group mean difference
#'
#'
#' @description
#' Computes a confidence interval for a population 2-group mean difference
#' using the sample means, sample standard deviations, and sample sizes. Use 
#' the t.test function for raw data input.
#'
#'  
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  m1     sample mean for group 1
#' @param  m2     sample mean for group 2
#' @param  sd1    sample standard deviation for group 1
#' @param  sd2    sample standard deviation for group 2
#' @param  n1     sample size for group 1
#' @param  n2     sample size for group 2
#'
#'
#' @return 
#' Returns a 2-row matrix. The columns are:
#' * Estimate - estimated mean difference
#' * SE - standard error
#' * t - t test statistic
#' * df - degrees of freedom
#' * p - p-value
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @references
#' \insertRef{Snedecor1980}{statpsych}
#'
#'
#' @examples
#' ci.mean2(.05, 15.4, 10.3, 2.67, 2.15, 30, 20)
#'
#' # Should return:
#' #                              Estimate       SE        t      df      
#' # Equal Variances Assumed:          5.1 1.602248 3.183029 48.0000 
#' # Equal Variances Not Assumed:      5.1 1.406801 3.625247 44.1137 
#' #                                          p       LL       UL
#' # Equal Variances Assumed:      0.0025578586 1.878465 8.321535
#' # Equal Variances Not Assumed:  0.0007438065 2.264986 7.935014
#'
#'
#' @importFrom stats qt
#' @importFrom stats pt
#' @export
ci.mean2 <- function(alpha, m1, m2, sd1, sd2, n1, n2) {
 df1 <- n1 + n2 - 2
 est <- m1 - m2
 v1 <- sd1^2
 v2 <- sd2^2
 vp <- ((n1 - 1)*v1 + (n2 - 1)*v2)/df1
 se1 <- sqrt(vp/n1 + vp/n2)
 t1 <- est/se1
 p1 <- 2*(1 - pt(abs(t1),df1))
 tcrit1 <- qt(1 - alpha/2, df1)
 ll1 <- est - tcrit1*se1
 ul1 <- est + tcrit1*se1
 se2 <- sqrt(v1/n1 + v2/n2)
 t2 <- est/se2
 df2 <- (se2^4)/(v1^2/(n1^3 - n1^2) + v2^2/(n2^3 - n2^2))
 p2 <- 2*(1 - pt(abs(t2),df2))
 tcrit2 <- qt(1 - alpha/2, df2)
 ll2 <- est - tcrit2*se2
 ul2 <- est + tcrit2*se2
 out1 <- t(c(est, se1, t1, df1, p1, ll1, ul1))
 out2 <- t(c(est, se2, t2, df2, p2, ll2, ul2))
 out <- rbind(out1, out2)
 colnames(out) <- c("Estimate", "SE", "t", "df", "p", "LL", "UL")
 rownames(out) <- c("Equal Variances Assumed:", "Equal Variances Not Assumed:")
 return(out)
}


#  ci.lc.mean.bs ==============================================================
#' Confidence interval for a linear contrast of means in a between-subjects
#' design
#' 
#'
#' @description
#' Computes a test statistic and confidence interval for a linear contrast
#' of means. This function computes both unequal variance (recommended) and 
#' equal variance confidence intervals and test statistics. A Satterthwaite 
#' adjustment to the degrees of freedom is used with the unequal variance 
#' method. 
#'
#'
#' @param     alpha  	alpha level for 1-alpha confidence
#' @param     m     	vector of group sample means
#' @param     sd    	vector of group sample standard deviations
#' @param     n     	vector of sample sizes
#' @param     v     	vector of betwen-subjects contrast coefficients
#' 
#'
#' @return 
#' Returns a 2-row matrix. The columns are:
#' * Estimate - estimated linear contrast
#' * SE - standard error
#' * t - t test statistic
#' * df - degrees of freedom
#' * p - p-value
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @references
#' \insertRef{Snedecor1980}{statpsych}
#'
#'
#' @examples
#' m <- c(33.5, 37.9, 38.0, 44.1)
#' sd <- c(3.84, 3.84, 3.65, 4.98)
#' n <- c(10,10,10,10)
#' v <- c(.5, .5, -.5, -.5)
#' ci.lc.mean.bs(.05, m, sd, n, v)
#'
#' # Should return:
#' #                              Estimate       SE         t       df 
#' # Equal Variances Assumed:        -5.35 1.300136 -4.114955 36.00000 
#' # Equal Variances Not Assumed:    -5.35 1.300136 -4.114955 33.52169 
#' #                                         p         LL        UL
#' # Equal Variances Assumed:     0.0002152581  -7.986797 -2.713203
#' # Equal Variances Not Assumed: 0.0002372436  -7.993583 -2.706417
#'
#'
#' @importFrom stats qt
#' @importFrom stats pt
#' @export
ci.lc.mean.bs <- function(alpha, m, sd, n, v) {
 est <- t(v)%*%m 
 k <- length(m)
 df1 <- sum(n) - k
 v1 <- sum((n - 1)*sd^2)/df1
 se1 <- sqrt(v1*t(v)%*%solve(diag(n))%*%v)
 t1 <- est/se1
 p1 <- 2*(1 - pt(abs(t1),df1))
 tcrit1 <- qt(1 - alpha/2, df1)
 ll1 <- est - tcrit1*se1
 ul1 <- est + tcrit1*se1
 v2 <- diag(sd^2)%*%(solve(diag(n)))
 se2 <- sqrt(t(v)%*%v2%*%v)
 t2 <- est/se2
 df2 <- (se2^4)/sum(((v^4)*(sd^4)/(n^2*(n - 1))))
 p2 <- 2*(1 - pt(abs(t2),df2))
 tcrit2 <- qt(1 - alpha/2, df2)
 ll2 <- est - tcrit2*se2
 ul2 <- est + tcrit2*se2
 out1 <- t(c(est, se1, t1, df1, p1, ll1, ul1))
 out2 <- t(c(est, se2, t2, df2, p2, ll2, ul2))
 out <- rbind(out1, out2)
 colnames(out) <- c("Estimate", "SE", "t", "df", "p", "LL", "UL")
 rownames(out) <- c("Equal Variances Assumed:", "Equal Variances Not Assumed:")
 return(out)
}


#  ci.tukey ===================================================================
#' Tukey-Kramer confidence intervals for all pairwise mean differences in a
#' betwen-subjects design
#' 
#'
#' @description
#' Computes heteroscedastic Tukey-Kramer (also known as Games-Howell) 
#' confidence intervals for all pairwise comparisons of population means 
#' using sample means, sample standard deviations, and samples sizes as 
#' input. A Satterthwaite adjustment to the degrees of freedom is used to 
#' improve the accuracy of the confidence intervals. 
#'
#'
#' @param  alpha   alpha level for simultaneous 1-alpha confidence
#' @param  m       vector of group sample means
#' @param  sd      vector of group sample standard deviations
#' @param  n       vector of sample sizes
#'
#'
#' @return 
#' Returns a matrix with the number of rows equal to the number
#' of pairwise comparisons. The columns are:
#' * Estimate - estimated mean difference
#' * SE - standard error
#' * t - t test statistic
#' * df - degrees of freedom
#' * p - p-value
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @references
#' \insertRef{Games1976}{statpsych}
#'
#'
#' @examples
#' m <- c(12.86, 17.57, 26.29)
#' sd <- c(3.185, 3.995, 3.773)
#' n <- c(20, 20, 20)
#' ci.tukey(.05, m, sd, n)
#'
#' # Should return:
#' #      Estimate       SE          t     p       df         LL         UL
#' # 1 2     -4.71 1.142459  -4.122686 0.001 36.20303  -7.501842  -1.918158
#' # 1 3    -13.43 1.104078 -12.163998 0.000 36.95915 -16.125713 -10.734287
#' # 2 3     -8.72 1.228730  -7.096758 0.000 37.87646 -11.717053  -5.722947
#'
#'
#' @importFrom stats qtukey
#' @importFrom stats ptukey
#' @importFrom utils combn
#' @export
ci.tukey <-function(alpha, m, sd, n) {
 a <- length(m)
 v1 <- sd^2/n
 v2 <- sd^4/(n^2*(n - 1))
 mean <- outer(m, m, '-')
 diff <- (-1)*mean[lower.tri(mean)]
 v1 <- outer(v1, v1, "+")
 v2 <- outer(v2, v2, "+")
 df = v1^2/v2
 df <- df[lower.tri(df)]
 SE <- sqrt(v1[lower.tri(v1)])
 t <- diff/SE
 q <- qtukey(p = 1 - alpha, nmeans = a, df = df)/sqrt(2)
 p <- 1 - ptukey(sqrt(2)*abs(t), nmeans = a, df = df)
 p <- round(p*1000)/1000
 LL <- diff - q*SE
 UL <- diff + q*SE
 pair <- t(combn(seq(1:a), 2))
 out <- cbind(pair, diff, SE, t, df, p, LL, UL)
 rownames(out) <- rep("", a*(a - 1)/2)
 return(out)
}


# ci.ratio.mean2 ==============================================================
#' Confidence interval for a 2-group mean ratio
#'
#'
#' Computes a confidence interval for a ratio of population means of 
#' ratio-scale measurements in a 2-group design. Equality of variances 
#' is not assumed.
#'
#'
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  y1     vector of scores for group 1
#' @param  y2     vector of scores for group 2
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Mean1 - estimated mean from group 1
#' * Mean2 - estimated mean from group 2
#' * Mean1/Mean2- estimated mean ratio
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @references
#' \insertRef{Bonett2020b}{statpsych}
#'
#'
#' @examples
#' y2 <- c(32, 39, 26, 35, 43, 27, 40, 37, 34, 29, 49, 42, 40)
#' y1 <- c(36, 44, 47, 42, 49, 39, 46, 31, 33, 48)
#' ci.ratio.mean2(.05, y1, y2)
#'
#' # Should return:
#' #
#' #      Mean1    Mean2 Mean1/Mean2        LL       UL
#' # [1,]  41.5 36.38462    1.140592 0.9897482 1.314425
#'
#'
#' @importFrom stats qt
#' @importFrom stats var
#' @export
ci.ratio.mean2 <- function(alpha, y1, y2){
 n1 <- length(y1)
 n2 <- length(y2)
 m1 <- mean(y1)
 m2 <- mean(y2)
 v1 <- var(y1)
 v2 <- var(y2)
 var <- v1/(n1*m1^2) + v2/(n2*m2^2)
 df <- var^2/(v1^2/(m1^4*(n1^3 - n1^2)) + v2^2/(m2^4*(n2^3 - n2^2)))
 tcrit <- qt(1 - alpha/2, df)
 est <- log(m1/m2)
 se <- sqrt(var)
 ll <- exp(est - tcrit*se)
 ul <- exp(est + tcrit*se)
 out <- t(c(m1, m2, exp(est), ll, ul))
 colnames(out) <- c("Mean1", "Mean2", "Mean1/Mean2", "LL", "UL")
 return(out)
}


#  ci.stdmean2 ================================================================
#' Confidence interval for a 2-group standardized mean difference
#' 
#'
#' @description
#' Computes confidence intervals for a population standardized mean difference. 
#' Unweighted, weighted, and single group variance standardizers are used. The 
#' square root weighted variance standardizer is recommended in 2-group 
#' nonexperimental designs with simple random sampling. The square root 
#' unweighted variance standardizer is recommended in 2-group experimental 
#' designs. The single group standard deviation standardizer can be used with 
#' experimental or nonexperimental designs. Equality of variances is not 
#' assumed.
#'
#'
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  m1     sample mean in group 1
#' @param  m2     sample mean in group 2
#' @param  sd1    sample standard deviation in group 1
#' @param  sd2    sample standard deviation in group 2
#' @param  n1     sample size in group 1
#' @param  n2     sample size in group 2
#'
#'
#' @return 
#' Returns a 4-row matrix. The columns are: 
#' * Estimate - bias adjusted standardized mean difference
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @references
#' \insertRef{Bonett2008}{statpsych}
#'
#'
#' @examples
#' ci.stdmean2(.05, 35.1, 26.7, 7.32, 6.98, 30, 30)
#' 
#' # Should return:
#' #                           Estimate        SE        LL       UL
#' # Unweighted standardizer:  1.159240 0.2844012 0.6170771 1.731909
#' # Weighted standardizer:    1.159240 0.2802826 0.6251494 1.723837
#' # Group 1 standardizer:     1.117605 0.2975582 0.5643375 1.730744
#' # Group 2 standardizer:     1.172044 0.3120525 0.5918268 1.815050
#'
#'
#' @importFrom stats qnorm
#' @export
ci.stdmean2 <- function(alpha, m1, m2, sd1, sd2, n1, n2) {
 z <- qnorm(1 - alpha/2)
 v1 <- sd1^2
 v2 <- sd2^2
 df1 <- n1 - 1
 df2 <- n2 - 1
 df3 <- n1 + n2 - 2
 adj1 <- 1 - 3/(4*df1 - 1)
 adj2 <- 1 - 3/(4*df2 - 1)
 adj3 <- 1 - 3/(4*df3 - 1)
 s <- sqrt((v1 + v2)/2)
 sp <- sqrt((df1*v1 + df2*v2)/df3)
 est1 <- (m1 - m2)/s
 est1u <- adj3*est1
 se1 <- sqrt(est1^2*(v1^2/df1 + v2^2/df2)/(8*s^4) + (v1/df1 + v2/df2)/s^2)
 ll1 <- est1 - z*se1
 ul1 <- est1 + z*se1
 est2 <- (m1 - m2)/sp
 est2u <- adj3*est2
 se2 <- sqrt(est2^2*(1/df1 + 1/df2)/8 + (v1/n1 + v2/n2)/sp^2)
 ll2 <- est2 - z*se2
 ul2 <- est2 + z*se2
 est3 <- (m1 - m2)/sd1
 est3u <- adj1*est3
 se3 <- sqrt(est3^2/(2*df1) + 1/df1 + v2/(df2*v1))
 ll3 <- est3 - z*se3
 ul3 <- est3 + z*se3
 est4 <- (m1 - m2)/sd2
 est4u <- adj2*est4
 se4 <- sqrt(est4^2/(2*df2) + 1/df2 + v1/(df1*v2))
 ll4 <- est4 - z*se4
 ul4 <- est4 + z*se4
 out1 <- t(c(est1u, se1, ll1, ul1))
 out2 <- t(c(est2u, se2, ll2, ul2))
 out3 <- t(c(est3u, se3, ll3, ul3))
 out4 <- t(c(est4u, se4, ll4, ul4))
 out <- rbind(out1, out2, out3, out4)
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 rownames1 <- c("Unweighted standardizer:", "Weighted standardizer:")
 rownames2 <- c("Group 1 standardizer:", "Group 2 standardizer:")
 rownames(out) <- c(rownames1, rownames2)
 return(out)
}
	

#  ci.stdmean.strat ===========================================================
#' Confidence interval for a 2-group standardized mean difference with 
#' stratified sampling
#'
#'
#' @description
#' Computes confidence intervals for a population standardized mean difference
#' in a 2-group nonexperimental design with stratified random sampling using a 
#' square root weighted variance standardizer or single group standard 
#' deviation standardizer.  Equality of variances is not assumed.
#'
#'
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  m1     sample mean in group 1
#' @param  m2     sample mean in group 2
#' @param  sd1    sample standard deviation in group 1
#' @param  sd2    sample standard deviation in group 2
#' @param  n1     sample size in group 1
#' @param  n2     sample size in group 2
#' @param  p1     proportion of total population in subpopulation 1
#'
#'
#' @return 
#' Returns a 3-row matrix. The columns are: 
#' * Estimate - bias adjusted standardized mean difference
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @references
#' \insertRef{Bonett2020a}{statpsych}
#'
#'
#' @examples
#' ci.stdmean.strat(.05, 30.2, 30.8, 10.5, 11.2, 200, 200, .533)
#'
#' # Should return:
#' #                           Estimate         SE         LL        UL
#' # Weighted standardizer: -0.05528428 0.10023259 -0.2518410 0.1410636
#' # Group 1 standardizer:  -0.05692722 0.10368609 -0.2603639 0.1460782
#' # Group 2 standardizer:  -0.05692722 0.09720571 -0.2440911 0.1369483
#'
#'
#' @importFrom stats qnorm
#' @export
ci.stdmean.strat <- function(alpha, m1, m2, sd1, sd2, n1, n2, p1) {
 z <- qnorm(1 - alpha/2)
 v1 <- sd1^2
 v2 <- sd2^2
 df1 <- n1 - 1
 df2 <- n2 - 1
 df3 <- n1 + n2 - 2
 adj1 <- 1 - 3/(4*df1 - 1)
 adj2 <- 1 - 3/(4*df2 - 1)
 adj3 <- 1 - 3/(4*df3 - 1)
 s <- sqrt(p1*v1 + (1 - p1)*v2)
 est1 <- (m1 - m2)/s
 est1u <- adj3*est1
 se1 <- sqrt(est1^2*(1/df1 + 1/df2)/8 + (v1/n1 + v2/n2)/s^2)
 ll1 <- est1 - z*se1
 ul1 <- est1 + z*se1
 est3 <- (m1 - m2)/sd1
 est3u <- adj1*est3
 se3 <- sqrt(est3^2/(2*df1) + 1/df1 + v2/(df2*v1))
 ll3 <- est3 - z*se3
 ul3 <- est3 + z*se3
 est4 <- (m1 - m2)/sd2
 est4u <- adj2*est3
 se4 <- sqrt(est4^2/(2*df2) + 1/df2 + v1/(df1*v2))
 ll4 <- est4 - z*se4
 ul4 <- est4 + z*se4
 out1 <- t(c(est1u, se1, ll1, ul1))
 out3 <- t(c(est3u, se3, ll3, ul3))
 out4 <- t(c(est4u, se4, ll4, ul4))
 out <- rbind(out1, out3, out4)
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 rownames1 <- c("Weighted standardizer:")
 rownames2 <- c("Group 1 standardizer:", "Group 2 standardizer:")
 rownames(out) <- c(rownames1, rownames2)
 return(out)
}


#  ci.lc.stdmean.bs ==========================================================
#' Confidence interval for a standardized linear contrast of means in a 
#' between-subjects design
#'
#'
#' @description
#' Computes confidence intervals for a population standardized linear contrast 
#' of means in a between-subjects design. The square root weighted variance 
#' standardizer is recommended in 2-group nonexperimental designs with simple 
#' random sampling. The square root unweighted variance standardizer is 
#' recommended in 2-group experimental designs. Equality of variances is not
#' assumed.
#'
#'
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  m      vector of group sample means
#' @param  sd     vector of group sample standard deviation
#' @param  n      vector of sample sizes
#' @param  v      vector of between-subjects contrast coefficients
#'
#'
#' @return 
#' Returns a 2-row matrix. The columns are:
#' * Estimate - bias adjusted standardized linear contrast
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @references
#' \insertRef{Bonett2008}{statpsych}
#'
#'
#' @examples
#' m <- c(33.5, 37.9, 38.0, 44.1)
#' sd <- c(3.84, 3.84, 3.65, 4.98)
#' n <- c(10,10,10,10)
#' v <- c(.5, .5, -.5, -.5)
#' ci.lc.stdmean.bs(.05, m, sd, n, v)
#'
#' # Should return:
#' #                           Estimate        SE        LL         UL
#' # Unweighted standardizer: -1.273964 0.3692800 -2.025039 -0.5774878
#' # Weighted standardizer:   -1.273964 0.3514511 -1.990095 -0.6124317
#'
#'
#' @importFrom stats qnorm
#' @export
ci.lc.stdmean.bs <- function(alpha, m, sd, n, v) {
 z <- qnorm(1 - alpha/2)
 var <- sd^2
 a <- length(m)
 s <- sqrt(sum(var)/a)
 df <- sum(n) - a
 adj <- 1 - 3/(4*df - 1)
 sp <- sqrt(sum((n - 1)*var)/df)
 est1 <- (t(v)%*%m)/s
 est2 <- (t(v)%*%m)/sp
 est1u <- adj*est1
 est2u <- adj*est2
 a1 <- est1^2/(2*a^2*s^4)
 a2 <- a1*sum((var^2/(n - 1)))
 a3 <- sum((v^2*var/(n - 1)))/s^2
 se1 <- sqrt(a2 + a3)
 ll1 <- est1 - z*se1
 ul1 <- est1 + z*se1
 a1 <- est2^2/(2*a^2)
 a2 <- a1*sum(1/(n - 1))
 a3 <- sum((v^2*var/n))/sp^2
 se2 <- sqrt(a2 + a3)
 ll2 <- est2 - z*se2
 ul2 <- est2 + z*se2
 out1 <- t(c(est1u, se1, ll1, ul1))
 out2 <- t(c(est2u, se2, ll2, ul2))
 out <- rbind(out1, out2)
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 rownames(out) <- c("Unweighted standardizer:", "Weighted standardizer:")
 return(out)
}


#  ci.mean.ps ================================================================
#' Confidence interval for a paired-samples mean difference
#'
#' @description
#' Computes a confidence interval for a population paired-samples mean 
#' difference using the sample means, sample standard deviations, sample 
#' correlation, and sample size. Use the t.test function for raw data input.
#'
#'  
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  m1     sample mean for measurement 1
#' @param  m2     sample mean for measurement 2
#' @param  sd1    sample standard deviation for measurement 1
#' @param  sd2    sample standard deviation for measurement 2
#' @param  cor    sample correlation between measurements
#' @param  n      sample size
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimated mean difference
#' * SE - standard error
#' * t - t test statistic
#' * df - degrees of freedom
#' * p - p-value
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' ci.mean.ps(.05, 58.2, 51.4, 7.43, 8.92, .537, 30)
#'
#' # Should return:
#' #      Estimate       SE        t df            p       LL       UL
#' # [1,]      6.8 1.455922 4.670578 29  6.33208e-05 3.822304 9.777696
#'
#'
#' @importFrom stats qt
#' @importFrom stats pt
#' @export
ci.mean.ps <- function(alpha, m1, m2, sd1, sd2, cor, n) {
 df <- n - 1
 tcrit <- qt(1 - alpha/2, df)
 vd <- sd1^2 + sd2^2 - 2*cor*sd1*sd2
 est <- m1 - m2
 se <- sqrt(vd/n)
 t <- est/se
 p <- 2*(1 - pt(abs(t), df))
 ll <- est - tcrit*se
 ul <- est + tcrit*se
 out <- t(c(est, se, t, df, p, ll, ul))
 colnames(out) <- c("Estimate", "SE", "t", "df", "p","LL", "UL")
 return(out)
}


#  ci.ratio.mean.ps ==========================================================
#' Confidence interval for a paired-samples mean ratio
#'
#'
#' @description
#' Compute a confidence interval for a ratio of population means of 
#' ratio-scale measurements in a paired-samples design. Equality of 
#' variances is not assumed.
#'
#'
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  y1     vector of measurement 1 scores
#' @param  y2     vector of measurement 2 scores
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Mean1 - estimated measurement 1 mean
#' * Mean2 - estimated measurement 2 mean
#' * Mean1/Mean2 - estimate of mean ratio
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @references
#' \insertRef{Bonett2020b}{statpsych}
#'
#'
#' @examples
#' y1 <- c(3.3, 3.6, 3.0, 3.1, 3.9, 4.2, 3.5, 3.3)
#' y2 <- c(3.0, 3.1, 2.7, 2.6, 3.2, 3.8, 3.2, 3.0)
#' ci.ratio.mean.ps(.05, y1, y2)
#'
#' # Should return:
#' #       Mean1 Mean2 Mean1/Mean2      LL       UL
#' # [1,] 3.4875 3.075    1.134146 1.09417 1.175583
#'
#'
#' @importFrom stats qt
#' @importFrom stats cor
#' @export
ci.ratio.mean.ps <- function(alpha, y1, y2){
 n <- length(y1)
 m1 <- mean(y1)
 m2 <- mean(y2)
 v1 <- var(y1)
 v2 <- var(y2)
 cor <- cor(y1,y2)
 var <- (v1/m1^2 + v2/m2^2 - 2*cor*sqrt(v1*v2)/(m1*m2))/n
 df <- n - 1
 tcrit <- qt(1 - alpha/2, df)
 est <- log(m1/m2)
 se <- sqrt(var)
 ll <- exp(est - tcrit*se)
 ul <- exp(est + tcrit*se)
 out <- t(c(m1, m2, exp(est), ll, ul))
 colnames(out) <- c("Mean1", "Mean2", "Mean1/Mean2", "LL", "UL")
 return(out)
}


#  ci.stdmean.ps =============================================================
#' Confidence interval for a paired-samples standardized mean difference
#'
#'
#' @description
#' Computes confidence intervals for a population standardized mean difference
#' in a paired-samples design. A square root unweighted variance standardizer 
#' and single measurement standard deviation standardizers are used. Equality 
#' of variances is not assumed.
#'
#'
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  m1     sample mean of measurement 1
#' @param  m2     sample mean of measurement 2
#' @param  sd1    sample standard deviation in measurement 1
#' @param  sd2    sample standard deviation in measurement 2
#' @param  cor    sample correlation between measurements
#' @param  n      sample size
#'
#'
#' @return 
#' Returns a 3-row matrix. The columns are:
#' * Estimate - bias adjusted standardized mean difference
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @references
#' \insertRef{Bonett2008}{statpsych}
#'
#'
#' @examples
#' ci.stdmean.ps(.05, 110.4, 102.1, 15.3, 14.6, .75, 25)
#'
#' # Should return:
#' #                              Estimate        SE        LL        UL
#' # Unweighted standardizer:     0.5433457 0.1609934 0.2394905 0.8705732
#' # Measurement 1 standardizer:  0.5253526 0.1615500 0.2258515 0.8591158
#' # Measurement 2 standardizer:  0.5505407 0.1692955 0.2366800 0.9003063
#'
#'
#' @importFrom stats qnorm
#' @export
ci.stdmean.ps <- function(alpha, m1, m2, sd1, sd2, cor, n) {
 z <- qnorm(1 - alpha/2)
 s <- sqrt((sd1^2 + sd2^2)/2)
 df <- n - 1
 v1 <- sd1^2
 v2 <- sd2^2
 adj1 <- sqrt((n - 2)/df)
 adj2 <- 1 - 3/(4*df - 1)
 vd <- v1 + v2 - 2*cor*sd1*sd2
 est1 <- (m1 - m2)/s
 est1u <- adj1*est1
 se1 <- sqrt(est1^2*(v1^2 + v2^2 + 2*cor^2*v1*v2)/(8*df*s^4) + vd/(df*s^2))
 ll1 <- est1 - z*se1
 ul1 <- est1 + z*se1
 est3 <- (m1 - m2)/sd1
 est3u <- adj2*est3
 se3 <- sqrt(est3^2/(2*df) + vd/(df*v1))
 ll3 <- est3 - z*se3
 ul3 <- est3 + z*se3
 est4 <- (m1 - m2)/sd2
 est4u <- adj2*est4
 se4 <- sqrt(est4^2/(2*df) + vd/(df*v2))
 ll4 <- est4 - z*se4
 ul4 <- est4 + z*se4
 out1 <- t(c(est1u, se1, ll1, ul1))
 out3 <- t(c(est3u, se3, ll3, ul3))
 out4 <- t(c(est4u, se4, ll4, ul4))
 out <- rbind(out1, out3, out4)
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 rownames1 <- c("Unweighted standardizer:")
 rownames2 <- c("Measurement 1 standardizer:", "Measurement 2 standardizer:")
 rownames(out) <- c(rownames1, rownames2)
 return(out)
}

#  ci.lc.stdmean.ws ==========================================================
#' Confidence interval for a standardized linear contrast of means in a
#' within-subjects design
#'
#'
#' @description
#' Computes a confidence interval for a population standardized linear 
#' contrast of means in a within-subjects design. A square root unweighted
#' variance standardizer is used. Equality of variances is not assumed but 
#' the correlations among the repeated measures are assumed to be 
#' approximately equal.
#'
#'
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  m      vector of within-subjects sample means
#' @param  sd     vector of within-subjects sample standard deviations
#' @param  cor    average correlation of measurement pairs
#' @param  n      sample size
#' @param  q      vector of within-subjects contrast coefficients
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Estimate - bias adjusted standardized linear contrast
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @references
#' \insertRef{Bonett2008}{statpsych}
#'
#'
#' @examples
#' m <- c(33.5, 37.9, 38.0, 44.1)
#' sd <- c(3.84, 3.84, 3.65, 4.98)
#' q <- c(.5, .5, -.5, -.5)
#' ci.lc.stdmean.ws(.05, m, sd, .672, 20, q)
#'
#' # Should return:
#' #       Estimate        SE        LL        UL
#' # [1,] -1.266557 0.1860897 -1.665992 -0.936534
#'
#'
#' @importFrom stats qnorm
#' @export
ci.lc.stdmean.ws <- function(alpha, m, sd, cor, n, q) {
 z <- qnorm(1 - alpha/2)
 a <- length(m)
 df <- n - 1
 adj <- sqrt((n - 2)/df)
 s <- sqrt(sum(sd^2)/a)
 est <- (t(q)%*%m)/s
 estu <- adj*est
 v1 <- est^2/(2*a^2*s^4*df)
 v2 <- sum(sd^4)
 v3 <- cor^2*t(sd^2)%*%sd^2 
 v4 <- sum(q^2*sd^2)
 v5 <- cor*t(q*sd)%*%(q*sd)
 se <- sqrt(v1*(v2 + v3) + (v4 - v5)/(df*s^2))
 ll <- est - z*se
 ul <- est + z*se
 out <- t(c(estu, se, ll, ul))
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 return(out)
}


#  ci.mad1 ====================================================================
#' Confidence interval for a single MAD
#'
#'
#' @description
#' Computes a confidence interval for a population mean absolute deviation 
#' from the median (MAD). 
#'
#'
#' @param  alpha   alpha level for 1-alpha confidence
#' @param  y       vector of scores
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * MAD - estimated mean absolute deviation
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @references
#' \insertRef{Bonett2003b}{statpsych}
#'
#'
#' @examples
#' y <- c(30, 20, 15, 10, 10, 60, 20, 25, 20, 30, 10, 5, 50, 40, 
#'        20, 10, 0, 20, 50)
#' ci.mad1(.05, y)
#'
#' # Should return:
#' #       MAD       LL       UL
#' # [1,] 12.5 7.962667 19.62282
#'
#'
#' @importFrom stats qnorm
#' @importFrom stats sd
#' @export
ci.mad1 <- function(alpha, y) {
 n <- length(y)
 z <- qnorm(1 - alpha/2)
 c <- n/(n - 1)
 median <- median(y)
 mad <- mean(abs(y - median))
 skw <- (mean(y) - median)/mad 
 kur <- (sd(y)/mad)^2
 se <- sqrt((skw^2 + kur - 1)/n)
 ll <- exp(log(c*mad) - z*se)
 ul <- exp(log(c*mad) + z*se)
 out <- t(c(c*mad, ll, ul))
 colnames(out) <- c("MAD", "LL", "UL")
 return(out)
}


#  ci.ratio.mad2 ==============================================================
#' Confidence interval for a 2-group MAD ratio
#'
#'
#' @description
#' Computes a confidence interval for a ratio of population MADs (mean absolute
#' deviation from median) in a 2-group design.
#'
#'
#' @param  alpha   alpha level for 1-alpha confidence
#' @param  y1      vector of scores for group 1
#' @param  y2      vector of scores for group 2
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * MAD1 - estimated MAD from group 1
#' * MAD2 - estimated MAD from group 2
#' * MAD1/MAD2 - estimate of MAD ratio
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#'
#' @references
#' \insertRef{Bonett2003b}{statpsych}
#'
#'
#' @examples
#' y1 <- c(32, 39, 26, 35, 43, 27, 40, 37, 34, 29)
#' y2 <- c(36, 44, 47, 42, 49, 39, 46, 31, 33, 48)
#' ci.ratio.mad2(.05, y1, y2)
#'
#' # Should return:
#' #          MAD1     MAD2  MAD1/MAD2        LL       UL
#' # [1,] 5.111111 5.888889  0.8679245 0.4520879 1.666253
#'
#'
#' @importFrom stats qnorm
#' @importFrom stats sd
#' @importFrom stats median
#' @export
ci.ratio.mad2 <- function(alpha, y1, y2) {
 z <- qnorm(1 - alpha/2)
 n1 <- length(y1)
 c1 <- n1/(n1 - 1)
 n2 <- length(y2)
 c2 <- n2/(n2 - 1)
 median1 <- median(y1)
 median2 <- median(y2)
 mad1 <- mean(abs(y1 - median1))
 skw1 <- (mean(y1) - median1)/mad1 
 kur1 <- (sd(y1)/mad1)^2
 var1 <- (skw1^2 + kur1 - 1)/n1
 mad2 <- mean(abs(y2 - median2))
 skw2 <- (mean(y2) - median2)/mad2 
 kur2 <- (sd(y2)/mad2)^2
 var2 <- (skw2^2 + kur2 - 1)/n2
 se <- sqrt(var1 + var2)
 c <- c1/c2
 est <- mad1/mad2
 ll <- exp(log(c*est) - z*se)
 ul <- exp(log(c*est) + z*se)
 out <- t(c(c1*mad1, c2*mad2, c*est, ll, ul))
 colnames(out) <- c("MAD1", "MAD2", "MAD1/MAD2", "LL", "UL")
 return(out)
}


#  ci.ratio.mad.ps ========================================================== 
#' Confidence interval for a paired-sample MAD ratio 
#'
#'
#' @description
#' Computes a confidence interval for a ratio of population MADs (mean absolute
#' deviation from median) in a paired-samples design.
#'
#' @param  alpha   alpha level for 1-alpha confidence
#' @param  y1      vector of measurement 1 scores
#' @param  y2      vector of measurement 2 scores
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * MAD1 - estimated MAD for measurement 1
#' * MAD2 - estimated MAD for measurement 2
#' * MAD1/MAD2 - estimate of MAD ratio
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#'
#' @references
#' \insertRef{Bonett2003a}{statpsych}
#'
#'
#' @examples
#' y2 <- c(21, 4, 9, 12, 35, 18, 10, 22, 24, 1, 6, 8, 13, 16, 19)
#' y1 <- c(67, 28, 30, 28, 52, 40, 25, 37, 44, 10, 14, 20, 28, 40, 51)
#' ci.ratio.mad.ps(.05, y1, y2)
#'
#' # Should return:
#' #          MAD1  MAD2  MAD1/MAD2       LL       UL
#' # [1,] 12.71429   7.5   1.695238 1.109176 2.590961
#'
#'
#' @importFrom stats qnorm
#' @importFrom stats sd
#' @importFrom stats median
#' @export
ci.ratio.mad.ps <- function(alpha, y1, y2) {
 z <- qnorm(1 - alpha/2)
 n <- length(y1); 
 c <- n/(n - 1)
 median1 <- median(y1);
 median2 <- median(y2);
 mad1 <- mean(abs(y1 - median1))
 skw1 <- (mean(y1) - median1)/mad1 
 kur1 <- (sd(y1)/mad1)^2
 var1 <- (skw1^2 + kur1 - 1)/n
 mad2 <- mean(abs(y2 - median2))
 skw2 <- (mean(y2) - median2)/mad2 
 kur2 <- (sd(y2)/mad2)^2
 var2 <- (skw2^2 + kur2 - 1)/n
 d1 <- abs(y1 - median1);
 d2 <- abs(y2 - median2)
 cor <- cor(d1, d2)
 se <- sqrt(var1 + var2 - 2*cor*sqrt(var1*var2))
 est <- mad1/mad2
 ll <- exp(log(est) - z*se)
 ul <- exp(log(est) + z*se)
 out <- t(c(c*mad1, c*mad2, est, ll, ul))
 colnames(out) <- c("MAD1", "MAD2", "MAD1/MAD2", "LL", "UL")
 return(out)
}


#  ci.cod1 =================================================================== 
#' Confidence interval for a single coefficient of dispersion
#'
#'
#' @description
#' Computes a confidence interval for a population coefficient of dispersion
#' (COD) which is defined as MAD/median. The COD is a robust alternative to 
#' the coefficient of variation and assumes ratio-scale scores.
#'
#' @param  alpha   alpha level for 1-alpha confidence
#' @param  y       vector of scores
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * COD - estimated coefficient of dispersion
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @references
#' \insertRef{Bonett2006}{statpsych}
#'
#'
#' @examples
#' y <- c(30, 20, 15, 10, 10, 60, 20, 25, 20, 30, 10, 5, 50, 40,
#'        20, 10, 0, 20, 50)
#' ci.cod1(.05, y)
#'
#' # Should return:
#' #            COD        LL       UL
#' # [1,] 0.5921053 0.3813259 1.092679
#'
#'
#' @importFrom stats qnorm
#' @importFrom stats var
#' @importFrom stats median
#' @export
ci.cod1 <-function(alpha, y) {
 z = qnorm(1 - alpha/2)
 n = length(y)
 c = n/(n - 1)
 a1 = round((n + 1)/2 - sqrt(n))
 b1 = n - a1 + 1
 med = median(y)
 m = mean(y)
 v = var(y)
 tau = mean(abs(y - med))
 del = (m - med)/tau
 gam = v/tau^2
 cod = tau/med
 y = sort(y)
 vareta = ((log(y[a1]) - log(y[b1]))/4)^2
 se1 = sqrt(vareta)
 vartau = (gam + del^2 - 1)/n
 se2 = sqrt(vartau)
 cov = del*se1/sqrt(n)
 k = sqrt(vareta + vartau - 2*cov)/(se1 + se2)
 a2 = round((n + 1)/2 - k*z*sqrt(n/4))
 b2 = n - a2 + 1
 L2 = log(y[a2])
 U2 = log(y[b2])
 L1 = log(c*tau) - k*z*se2
 U1 = log(c*tau) + k*z*se2
 LL = exp(L1 - U2)
 UL = exp(U1 - L2)
 out = t(c(cod, LL, UL))
 colnames(out) = c("COD", "LL", "UL")
 return(out)
}


#  ci.median1 =================================================================
#' Confidence interval for a single median 
#'
#'
#' @description
#' Computes a confidence interval for a single population median.
#'
#'
#' @param  alpha   alpha level for 1-alpha confidence
#' @param  y       vector of scores
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Median - estimated median
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @references
#' \insertRef{Snedecor1980}{statpsych}
#'
#'
#' @examples
#' y <- c(30, 20, 15, 10, 10, 60, 20, 25, 20, 30, 10, 5, 50, 40,
#'        20, 10, 0, 20, 50)
#' ci.median1(.05, y)
#'
#' # Should return:
#' #      Median       SE LL UL
#' # [1,]     20 5.390263 10 30
#'
#'
#' @importFrom stats qnorm
#' @importFrom stats pbinom
#' @importFrom stats median
#' @export
ci.median1 <- function(alpha, y) {
 n <- length(y)
 y <- sort(y)
 z <- qnorm(1 - alpha/2)
 median <- median(y)
 c1 <- round((n - z*sqrt(n))/2)
 if (c1 < 1) {c1 = 1}
 ll <- y[c1]
 ul <- y[n - c1 + 1]
 a <- round((n + 1)/2 - sqrt(n))
 if (a < 1) {a = 1}
 ll1 <- y[a]
 ul1 <- y[n - a + 1]
 p <- pbinom(a - 1, size = n, prob = .5)
 z0 <- qnorm(1 - p)
 se <- (ul1 - ll1)/(2*z0)
 out <- t(c(median, se, ll, ul))
 colnames(out) <- c("Median", "SE", "LL", "UL")
 return(out)
}


#  ci.median2 =================================================================
#' Confidence interval for a 2-group median difference
#'
#'
#' @description
#' Computes a confidence interval for a 2-group difference of population 
#' medians.
#'
#'
#' @param  alpha   alpha level for 1-alpha confidence
#' @param  y1      vector of scores for group 1
#' @param  y2      vector of scores for group 2
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Median1 - estimated median from group 1
#' * Median2 - estimated median from group 2
#' * Median1-Median2 - estimated difference in medians
#' * SE - standard error of the difference
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @references
#' \insertRef{Bonett2002}{statpsych}
#'
#'
#' @examples
#' y1 <- c(32, 39, 26, 35, 43, 27, 40, 37, 34, 29)
#' y2 <- c(36, 44, 47, 42, 49, 39, 46, 31, 33, 48)
#' ci.median2(.05, y1, y2)
#'
#' # Should return:
#' #       Median1 Median2 Median1-Median2       SE        LL          UL
#' # [1,]     34.5      43            -8.5 4.316291 -16.95977 -0.04022524
#'
#'
#' @importFrom stats qnorm
#' @importFrom stats median
#' @importFrom stats pbinom
#' @export
ci.median2 <- function(alpha, y1, y2) {
 z <- qnorm(1 - alpha/2)
 n1 <- length(y1)
 y1 <- sort(y1)
 n2 <- length(y2)
 y2 <- sort(y2)
 median1 <- median(y1)
 median2 <- median(y2)
 a1 <- round(n1/2 - sqrt(n1))
 if (a1 < 1) {a1 = 1}
 l1 <- y1[a1]
 u1 <- y1[n1 - a1 + 1]
 p <- pbinom(a1 - 1, size = n1, prob = .5)
 z0 <- qnorm(1 - p)
 se1 <- (u1 - l1)/(2*z0)
 a2 <- round(n2/2 - sqrt(n2))
 if (a2 < 1) {a2 = 1}
 l2 <- y2[a2]
 u2 <- y2[n2 - a2 + 1]
 p <- pbinom(a2 - 1, size = n2, prob = .5)
 z0 <- qnorm(1 - p)
 se2 <- (u2 - l2)/(2*z0)
 diff <- median1 - median2
 se <- sqrt(se1^2 + se2^2)
 ll <- diff - z*se
 ul <- diff + z*se
 out <- t(c(median1, median2, diff, se, ll, ul))
 colnames(out) <- c("Median1", "Median2", "Median1-Median2", "SE", "LL", "UL")
 return(out)
}


#  ci.ratio.median2 ==========================================================
#' Confidence interval for a 2-group median ratio
#'
#'
#' @description
#' Computes a confidence interval for a ratio of population medians of 
#' ratio-scale measurements in a 2-group design.
#' 
#'
#' @param  alpha   alpha level for 1-alpha confidence
#' @param  y1      vector of scores for group 1
#' @param  y2      vector of scores for group 2
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Median1 - estimated median from group 1
#' * Median2 - estimated median from group 2
#' * Median1/Median2 - estimated ratio of medians
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @references
#' \insertRef{Bonett2020b}{statpsych}
#'
#'
#' @examples
#' y2 <- c(32, 39, 26, 35, 43, 27, 40, 37, 34, 29, 49, 42, 40)
#' y1 <- c(36, 44, 47, 42, 49, 39, 46, 31, 33, 48)
#' ci.ratio.median2(.05, y1, y2)
#'
#' # Should return:
#' #      Median1 Median2 Median1/Median2       LL       UL
#' # [1,]      43      37        1.162162 0.927667 1.455933
#'
#'
#' @importFrom stats qnorm
#' @importFrom stats pbinom
#' @importFrom stats median
#' @export
ci.ratio.median2 <- function(alpha, y1, y2) {
 z <- qnorm(1 - alpha/2)
 n1 <- length(y1)
 y1 <- sort(y1)
 n2 <- length(y2)
 y2 <- sort(y2)
 median1 <- median(y1)
 median2 <- median(y2)
 a1 <- round(n1/2 - sqrt(n1))
 if (a1 < 1) {a1 = 1}
 l1 <- log(y1[a1])
 u1 <- log((y1[n1 - a1 + 1]))
 p <- pbinom(a1 - 1, size = n1, prob = .5)
 z0 <- qnorm(1 - p)
 se1 <- (u1 - l1)/(2*z0)
 a2 <- round(n2/2 - sqrt(n2))
 if (a2 < 1) {a2 = 1}
 l2 <- log(y2[a2])
 u2 <- log((y2[n2 - a2 + 1]))
 p <- pbinom(a2 - 1, size = n2, prob = .5)
 z0 <- qnorm(1 - p)
 se2 <- (u2 - l2)/(2*z0)
 se <- sqrt(se1^2 + se2^2)
 diff <- log(median1) - log(median2)
 ll <- exp(diff - z*se)
 ul <- exp(diff + z*se)
 out <- t(c(median1, median2, exp(diff), ll, ul))
 colnames(out) <- c("Median1", "Median2", "Median1/Median2", "LL", "UL")
 return(out)
}


#  ci.lc.median.bs ===========================================================
#' Confidence interval for a linear contrast of medians in a between-subjects 
#' design
#'
#'
#' @description
#' Computes a confidence interval for a linear contrast of medians in a  
#' between-subjects design using sample medians and their standard errors.
#' The sample median and standard error for each group can be computed using
#' the ci.median1 function.
#'
#'
#' @param  alpha   alpha level for 1-alpha confidence
#' @param  m       vector of group sample medians
#' @param  se      vector of group standard errors 
#' @param  v       vector of between-subjects contrast coefficients
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimated linear contrast of medians
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @references
#' \insertRef{Bonett2002}{statpsych}
#'
#'
#' @examples
#' m <- c(46.13, 29.19, 30.32, 49.15)
#' se <- c(6.361, 5.892, 4.887, 6.103)
#' v <- c(1, -1, -1, 1)
#' ci.lc.median.bs(.05, m, se, v)
#'
#' # Should return:
#' #      Estimate       SE       LL       UL
#' # [1,]    35.77 11.67507 12.88727 58.65273
#'
#'
#' @importFrom stats qnorm
#' @export
ci.lc.median.bs <- function(alpha, m, se, v) {
 est <- t(v)%*%m
 se <- sqrt(t(v)%*%diag(se^2)%*%v)
 zcrit <- qnorm(1 - alpha/2)
 ll <- est - zcrit*se
 ul <- est + zcrit*se
 out <- t(c(est, se, ll, ul))
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 return(out)
}


#  ci.median.ps ==============================================================
#' Confidence interval for a paired-samples median difference
#'
#'
#' @description
#' Computes a confidence interval for a difference of population medians in 
#' a paired-samples design. This function also computes the standard errors
#' for each median and the covariance between the two medians.
#'
#'
#' @param  alpha   alpha level for 1-alpha confidence
#' @param  y1      vector of scores for measurement 1
#' @param  y2      vector of scores for measurement 2
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Median1 - estimated median for measurement 1
#' * Median2 - estimated median for measurement 2
#' * Median1-Median2 - estimated difference of medians
#' * SE1 - standard error of median 1
#' * SE2 - standard error of median 2
#' * cov - covariance of medians
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @references
#' \insertRef{Bonett2020}{statpsych}
#'
#'
#' @examples
#' y1 <- c(21, 4, 9, 12, 35, 18, 10, 22, 24, 1, 6, 8, 13, 16, 19)
#' y2 <- c(67, 28, 30, 28, 52, 40, 25, 37, 44, 10, 14, 20, 28, 40, 51)
#' ci.median.ps(.05, y1, y2)
#'
#' # Should return:
#' #    Median1 Median2 Median1-Median2       SE        LL        UL  
#' # [1,]    13      30             -17 3.362289 -23.58996 -10.41004 
#' #        SE1      SE2      cov
#' #   3.085608 4.509735 9.276849
#'
#'
#' @importFrom stats qnorm
#' @importFrom stats pbinom
#' @importFrom stats median
#' @export
ci.median.ps <- function(alpha, y1, y2) {
 z <- qnorm(1 - alpha/2)
 n <- length(y1)
 median1 <- median(y1)
 median2 <- median(y2)
 a1 <- (y1 < median1)
 a2 <- (y2 < median2)
 a3 <- a1 + a2
 a4 <- sum(a3 == 2)
 a <- round(n/2 - sqrt(n))
 if (a < 1) {a = 1}
 p <- pbinom(a - 1, size = n, prob = .5)
 z0 <- qnorm(1 - p)
 y1 <- sort(y1)
 y2 <- sort(y2)
 L1 <- y1[a]
 U1 <- y1[n - a + 1]
 se1 <- (U1 - L1)/(2*z0)
 L2 <- y2[a]
 U2 <- y2[n - a + 1]
 se2 <- (U2 - L2)/(2*z0)
 if (n/2 == trunc(n/2)) {
   p00 <- (sum(a4) + .25)/(n + 1)
 } else {
   p00 <- (sum(a4) + .25)/n 
 }
 cov <- (4*p00 - 1)*se1*se2
 diff <- median1 - median2
 se <- sqrt(se1^2 + se2^2 - 2*cov)
 ll <- diff - z*se
 ul <- diff + z*se
 out <- t(c(median1, median2, diff, se, ll, ul, se1, se2, cov))
 colnames(out) <- c("Median1", "Median2", "Median1-Median2", "SE", "LL", "UL", "SE1", "SE2", "cov")
 return(out)
}


#  ci.ratio.median.ps ========================================================
#' Confidence interval for a paired-samples median ratio
#'
#'
#' @description
#' Computes a confidence interval for a ratio of population medians in a 
#' paired-samples design.
#'
#'
#' @param  alpha   alpha level for 1-alpha confidence
#' @param  y1      vector of scores for measurement 1
#' @param  y2      vector of scores for measurement 2
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Median1 - estimated median from measurement 1
#' * Median2 - estimated median from measurement 2
#' * Median1/Median2 - estimated ratio of medians
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @references
#' \insertRef{Bonett2020b}{statpsych}
#'
#'
#' @examples
#' y1 <- c(21, 4, 9, 12, 35, 18, 10, 22, 24, 1, 6, 8, 13, 16, 19)
#' y2 <- c(67, 28, 30, 28, 52, 40, 25, 37, 44, 10, 14, 20, 28, 40, 51)
#' ci.ratio.median.ps(.05, y1, y2)
#'
#' # Should return:
#' #         Median1  Median2   Median1/Median2        LL        UL
#' # [1,]         13       30         0.4333333 0.3094838 0.6067451
#'
#'
#' @importFrom stats qnorm
#' @importFrom stats pbinom
#' @importFrom stats median
#' @export
ci.ratio.median.ps <- function(alpha, y1, y2) {
 z <- qnorm(1 - alpha/2)
 n <- length(y1)
 med1 <- median(y1)
 med2 <- median(y2)
 a1 <- (y1 < med1)
 a2 <- (y2 < med2)
 a3 <- a1 + a2
 f00 <- sum(a3 == 2)
 if (n/2 == trunc(n/2)) {
   p00 <- (f00 + .25)/(n + 1)
 } else {
   p00 <- (f00 + .25)/n 
 }
 o1 <- round(n/2 - sqrt(n))
 if (o1 < 1) {o1 = 1}
 o2 <- n - o1 + 1
 p <- pbinom(o1 - 1, size = n, prob = .5)
 z0 <- qnorm(1 - p)
 y1 <- sort(y1)
 y2 <- sort(y2)
 l1 <- log(y1[o1])
 u1 <- log(y1[o2])
 se1 <- (u1 - l1)/(2*z0)
 l2 <- log(y2[o1])
 u2 <- log(y2[o2])
 se2 <- (u2 - l2)/(2*z0)
 cov <- (4*p00 - 1)*se1*se2
 logratio <- log(med1/med2)
 se <- sqrt(se1^2 + se2^2 - 2*cov)
 ll <- exp(logratio - z*se)
 ul <- exp(logratio + z*se)
 out <- t(c(med1, med2, exp(logratio), ll, ul))
 colnames(out) <- c("Median1", "Median2", "Median1/Median2", "LL", "UL")
 return(out)
}


#  ci.mann ====================================================================
#' Confidence interval for a Mann-Whitney probability
#'
#'
#' @description
#' Computes a distribution-free confidence interval for the Mann-Whitney 
#' parameter. In a 2-group experiment, this parameter is the proportion of 
#' members in the population with scores that would be higher under treatment 1
#' than treatment 2. In a 2-group nonexperiment where participants are sampled 
#' from two subpopulations of sizes N1 and N2, the parameter is the proportion 
#' of all N1 x N2 pairs in which a member from subpopulation 1 has a larger
#' score than a member from subpopulation 2.
#'
#'
#' @param  alpha   alpha level for 1-alpha confidence
#' @param  y1      vector of scores for group 1
#' @param  y2      vector of scores for group 2
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimate of probability
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @references
#' \insertRef{Sen1967}{statpsych}
#'
#'
#' @examples
#' y1 <- c(32, 39, 26, 35, 43, 27, 40, 37, 34, 29)
#' y2 <- c(36, 44, 47, 42, 49, 39, 46, 31, 33, 48)
#' ci.mann(.05, y1, y2)
#'
#' # Should return:
#' #      Estimate        SE        LL UL
#' # [1,]    0.795 0.1401834 0.5202456  1
#'
#'
#' @importFrom stats qnorm
#' @export
ci.mann <- function(alpha, y1, y2){
 z <- qnorm(1 - .05/2)
 y <- c(y1,y2)
 n1 <- length(y1)
 n2 <- length(y2)
 n <- n1 + n2
 r <- rank(y)
 r1 <- r[1:n1]
 r2 <- r[(n1 + 1):n]
 m1 <- mean(r1)
 m2 <- mean(r2)
 seq1 <- seq(1,n1,1)
 seq2 <- seq(1,n2,1)
 a1 <- sum((r1 - seq1)^2)
 a2 <- sum((r2 - seq2)^2)
 v1 <- (a1 - n1*(m1 - (n1 + 1)/2)^2)/((n1 - 1)*n^2)
 v2 <- (a2 - n2*(m2 - (n2 + 1)/2)^2)/((n2 - 1)*n^2)
 u <- sum(r2) - n2*(n2 + 1)/2
 est <- u/(n1*n2)
 se <- sqrt((n2*v1 + n1*v2)/(n1*n2))
 ll <- est - z*se
 ul <- est + z*se
 if (ul > 1) {ul = 1}
 if (ll < 0) {ll = 0}
 out <- t(c(est, se, ll, ul))
 colnames(out) <- c("Estimate", "SE",  "LL", "UL")	
 return(out)
}


#  ci.random.anova1 ===========================================================
#' Confidence intervals for parameters of one-way random effects ANOVA
#'
#'
#' @description
#' Computes estimates and confidence intervals for four parameters of the
#' one-way random effects ANOVA: 1) the superpopulation grand mean, 2) the 
#' square-root within-group variance component, 3) the square-root 
#' between-group variance component, and 4) the omega-squared coefficient. 
#' This function assumes equal sample sizes.
#'
#'
#' @param   alpha  1 - alpha confidence 
#' @param   m      vector of group sample means 
#' @param   sd     vector of group standard deviations 
#' @param   n      sample size per group
#'
#'
#' @return 
#' Returns a 4-row matrix. The rows are:
#' * Grand mean - the mean of the superpopulation of means
#' * Within SD - the square-root within-group variance component
#' * Between SD - the square-root between-group variance component
#' * Omega-squared - the omega-squared coefficient
#'
#'
#' The columns are:
#' * Estimate - estimate of parameter
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @examples
#' m <- c(56.1, 51.2, 60.3, 68.2, 48.9, 70.5)
#' sd <- c(9.45, 8.79, 9.71, 8.90, 8.31, 9.75)
#' ci.random.anova1(.05, m, sd, 20)
#'
#' # Should return:
#' #                 Estimate         LL         UL
#' # Grand mean     59.200000 49.9363896 68.4636104
#' # Within SD:      9.166782  8.0509046 10.4373219
#' # Between SD:     8.585948  8.3239359  8.8562078
#' # Omega-squared:  0.467317  0.2284142  0.8480383
#'
#'
#' @importFrom stats qf
#' @importFrom stats qnorm
#' @importFrom stats qt
#' @export
ci.random.anova1 <- function(alpha, m, sd, n) {
 z <- qnorm(1 - alpha/2)
 a <- length(m)
 t <- qt(1 - alpha/2, a - 1)
 nt <- n*a
 dfe <- nt - a
 dfb <- a - 1
 v <- sd^2
 grandmean <- sum(n*m)/nt
 SSe <- sum((n - 1)*v)
 MSe <- SSe/dfe
 SSb <- sum(n*(m - grandmean)^2)
 MSb <- SSb/dfb
 F <- MSb/MSe
 varb <- (MSb - MSe)/n
 osqr <- (MSb - MSe)/(MSb + (n - 1)*MSe)
 se.mean <- sqrt(MSb/(a*n))
 se.vare <- sqrt(2*MSe^2/((dfb + 1)*(n - 1)))
 se.varb <- sqrt(2*(MSe^2/((dfb + 1)*(n - 1)) + MSb^2/dfb)/n^2)
 LL.m <- grandmean - t*se.mean
 UL.m <- grandmean + t*se.mean
 LL.e <- sqrt(exp(log(MSe) - z*se.vare/MSe))
 UL.e <- sqrt(exp(log(MSe) + z*se.vare/MSe))
 LL.b <- sqrt(exp(log(varb) - z*se.varb/MSb))
 UL.b <- sqrt(exp(log(varb) + z*se.varb/MSb))
 F1 <- qf(1 - alpha/2, dfb, dfe)
 F2 <- qf(alpha/2, dfb, dfe)
 LL.o <- (F/F1 - 1)/(n + F/F1 - 1)
 UL.o <- (F/F2 - 1)/(n + F/F2 - 1)
 out1 <- t(c(grandmean, LL.m, UL.m))
 out2 <- t(c(sqrt(MSe), LL.e, UL.e))
 out3 <- t(c(sqrt(varb), LL.b, UL.b))
 out4 <- t(c(osqr, LL.o, UL.o))
 out <- rbind(out1, out2, out3, out4)
 rownames(out) <- c("Grand mean", "Within SD:", "Between SD:", "Omega-squared:")
 colnames(out) = c("Estimate", "LL", "UL")
 return(out)
}


#  ci.cronbach ===============================================================
#' Confidence interval for a Cronbach reliability
#'
#'
#' @description
#' Computes a confidence interval for a population Cronbach reliability. 
#' The point estimate of Cronbach's reliability assumes essentially 
#' tau-equivalent measurements and this confidence interval assumes parallel
#' measurements. 
#'
#'  
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  rel    sample value of Cronbach's reliability  
#' @param  a      number of measurements (items, raters, etc.)
#' @param  n	  sample size
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @references
#' \insertRef{Feldt1965}{statpsych}
#'
#'
#' @examples
#' ci.cronbach(.05, .85, 7, 89)
#'
#' # Should return:
#' #              LL        UL
#' # [1,]  0.7971254 0.8931436   
#'  
#' 
#' @importFrom stats qf
#' @export
ci.cronbach <- function(alpha, rel, a, n) {
 df1 = n - 1
 df2 = n*(a - 1)
 f1 = qf(1 - alpha/2, df1, df2)
 f2 = qf(1 - alpha/2, df2, df1)
 f0 <- 1/(1 - rel)
 ll <- 1 - f1/f0
 ul <- 1 - 1/(f0*f2)
 out <- t(c(ll, ul))
 colnames(out) = c("LL", "UL")
 return(out)
}


# ================== File 1: Sample Size for Desired Precision ================
# size.ci.mean1 ===============================================================
#' Sample size for a single mean confidence interval
#'
#'
#' @description
#' Computes the sample size required to estimate a single population mean with
#' desired confidence interval precision. Set the variance planning value to   
#' the largest value within a plausible range for a conservatively large  
#' sample size. 
#'
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  var    planning value of response variable variance
#' @param  w      desired confidence interval width
#'
#'
#' @return 
#' Returns the required sample size
#'
#'
#' @examples
#' size.ci.mean1(.05, 264.4, 10)
#'
#' # Should return:
#' #      Sample size
#' # [1,]          43
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.ci.mean1 <- function(alpha, var, w) {
 z <- qnorm(1 - alpha/2)
 n <- ceiling(4*var*(z/w)^2 + z^2/2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 return(out)
}


#  size.ci.mean2  ============================================================ 
#' Sample size for a 2-group mean difference confidence interval
#'
#'
#' @description
#' Computes the sample size for each group required to estimate a population 
#' mean difference with desired confidence interval precision in a 2-group 
#' design. Set R = 1 for equal sample sizes. Set the variance planning value   
#' to the largest value within a plausible range for a conservatively large  
#' sample size. 
#'
#'
#' @param  alpha  alpha level for 1-alpha confidence 
#' @param  var    planning value of average within-group variance
#' @param  w      desired confidence interval width
#' @param  R      n2/n1 ratio 
#'
#'
#' @return 
#' Returns the required sample size for each group
#'
#'
#' @examples
#' size.ci.mean2(.05, 37.1, 5, 1)
#'
#' # Should return:
#' #      n1  n2
#' # [1,] 47  47
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.ci.mean2 <- function(alpha, var, w, R) {
 z <- qnorm(1 - alpha/2)
 n1 <- ceiling(4*var*(1 + 1/R)*(z/w)^2 + z^2/4)
 n2 <- ceiling(R*n1)
 out <- t(c(n1, n2))
 colnames(out) <- c("n1", "n2")
 return(out)
}


#  size.ci.stdmean2 ========================================================== 
#' Sample size for a 2-group standardized mean difference confidence interval
#'
#'
#' @description
#' Computes the sample size in each group required to estimate a population 
#' standardized mean difference with desired confidence interval precision 
#' in a 2-group design. Set the standardized mean difference planning value 
#' to the largest value within a plausible range for a conservatively large 
#' sample size. Set R = 1 for equal sample sizes.
#'
#'
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  d      planning value of standardized mean difference
#' @param  w      desired confidence interval width
#' @param  R      n2/n1 ratio
#'
#'
#' @references
#' \insertRef{Bonett2009}{statpsych}
#'
#'
#' @return 
#' Returns the required sample size for each group
#'
#'
#' @examples
#' size.ci.stdmean2(.05, .75, .5, 1)
#'
#' # Should return:
#' #       n1   n2
#' # [1,] 132  132
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.ci.stdmean2 <- function(alpha, d, w, R) {
 z <- qnorm(1 - alpha/2)
 n1 <- ceiling((d^2*(1 + R)/(2*R) + 4*(1 + R)/R)*(z/w)^2)
 n2 <- ceiling(R*n1)
 out <- t(c(n1, n2))
 colnames(out) <- c("n1", "n2")
 return(out)
}


#  size.ci.ratio.mean2 ======================================================= 
#' Sample size for a 2-group mean ratio confidence interval
#'
#'
#' @description
#' Computes the sample size in each group required to estimate a ratio of 
#' population means with desired confidence interval precision in a 2-group
#' design. This function requires planning values for each mean and the sample 
#' size requirement is very sensitive to these planning values. Set the 
#' variance planning value to the largest value within a plausible range for a
#' conservatively large sample size. 
#'
#'
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  var    planning value of average within-group variance
#' @param  m1     planning value of mean for group 1
#' @param  m2     planning value of mean for group 2
#' @param  r      desired upper to lower confidence interval endpoint ratio
#' @param  R      n2/n1 ratio
#'
#'
#' @return 
#' Returns the required sample size for each group
#'
#'
#' @examples
#' size.ci.ratio.mean2(.05, .4, 3.5, 3.1, 1.2, 2)
#'
#' # Should return:
#' #      n1   n2
#' # [1,] 53  106
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.ci.ratio.mean2 <- function(alpha, var, m1, m2, r, R) {
 z <- qnorm(1 - alpha/2)
 n1 <- ceiling(4*var*(1 + 1/R)*(1/m1^2 + 1/m2^2)*(z/log(r))^2 + z^2/4)
 n2 <- ceiling(R*n1)
 out <- t(c(n1, n2))
 colnames(out) <- c("n1", "n2")
 return(out)
}


#  size.ci.lc.mean.bs ========================================================
#' Sample size for a between-subjects mean linear contrast confidence interval
#'
#'
#' @description
#' Computes the sample size in each group (assuming equal sample sizes) 
#' required to estimate a linear contrast of population means with desired 
#' confidence interval precision in a between-subjects design. Set the 
#' variance planning value to the largest value within a plausible range
#' for a conservatively large sample size. 
#'
#'
#' @param  alpha  alpha level for 1-alpha confidence 
#' @param  var    planning value of average within-group variance  
#' @param  w      desired confidence interval width
#' @param  v      vector of between-subjects contrast coefficients 
#'
#'
#' @return 
#' Returns the required sample size for each group
#'
#'
#' @examples
#' v <- c(.5, .5, -1)
#' size.ci.lc.mean.bs(.05, 5.62, 2.0, v)
#'
#' # Should return:
#' #      Sample size per group
#' # [1,]                    34
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.ci.lc.mean.bs <- function(alpha, var, w, v) {
 z <- qnorm(1 - alpha/2)
 m <- length(v) - sum(v == 0)
 n <- ceiling(4*var*(t(v)%*%v)*(z/w)^2 + z^2/(2*m))
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size per group"
 return(out)
}


#  size.ci.lc.stdmean.bs ====================================================== 
#' Sample size for a between-subjects standardized mean linear contrast 
#' confidence interval
#'
#'
#' @description
#' Computes the sample size in each group (assuming equal sample sizes) 
#' required to estimate a standardized linear contrast of population means 
#' with desired confidence interval precision in a between-subjects design.
#' Set the standardized mean contrast planning value to the largest value 
#' within a plausible range for a conservatively large sample size. 
#'
#'
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  d      planning value of standardized linear contrast  
#' @param  w      desired confidence interval width
#' @param  v      vector of between-subjects contrast coefficients 
#'
#'
#' @references
#' \insertRef{Bonett2009}{statpsych}
#'
#'
#' @return 
#' Returns the required sample size for each group
#'
#'
#' @examples
#' v <- c(.5, .5, -.5, -.5)
#' size.ci.lc.stdmean.bs(.05, 1, .6, v)
#'
#' # Should return:
#' #      Sample size per group
#' # [1,]                    49
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.ci.lc.stdmean.bs <- function(alpha, d, w, v) {
 z <- qnorm(1 - alpha/2)
 a <- length(v)
 n <- ceiling((2*d^2/a + 4*(t(v)%*%v))*(z/w)^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size per group"
 return(out)	
}


#  size.ci.mean.ps ============================================================
#' Sample size for a paired-samples mean difference confidence interval
#'
#'
#' @description
#' Computes the sample size required to estimate a difference in population  
#' means with desired confidence interval precision in a paired-samples 
#' design. This function requires a planning value for the average of the  
#' variances for the two measurements. Set the correlation planning value to   
#' the smallest value within a plausible range for a conservatively large  
#' sample size. Set the variance planning value to the largest value within
#' a plausible range for a conservatively large sample size. 
#'
#'
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  var    planning value of average variance of the two measurements
#' @param  cor    planning value of correlation between measurements
#' @param  w      desired confidence interval width
#'
#'
#' @return 
#' Returns the required sample size
#'
#'
#' @examples
#' size.ci.mean.ps(.05, 265, .8, 10)
#'
#' # Should return:
#' #      Sample size
#' # [1,]          19
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.ci.mean.ps <- function(alpha, var, cor, w) {
 z <- qnorm(1 - alpha/2)
 n <- ceiling(8*(1 - cor)*var*(z/w)^2 + z^2/2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 return(out)	
}


#  size.ci.stdmean.ps ========================================================
#' Sample size for a paired-samples standardized mean difference confidence 
#' interval
#'
#'
#' @description
#' Computes the sample size required to estimate a population standardized 
#' mean difference with desired confidence interval precision in a 
#' paired-samples design. Set the standardized mean difference planning value 
#' to the largest value within a plausible range for a conservatively large 
#' sample size. Set the correlation planning value to the smallest value
#' within a plausible range for a conservatively large sample size.
#'
#'
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  d      planning value of standardized mean difference  
#' @param  cor    planning value of correlation between measurements
#' @param  w      desired confidence interval width
#'
#'
#' @references
#' \insertRef{Bonett2009}{statpsych}0
#'
#'
#' @return 
#' Returns the required sample size 
#'
#'
#' @examples
#' size.ci.stdmean.ps(.05, 1, .65, .6)
#'
#' # Should return:
#' #      Sample size
#' # [1,]          46
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.ci.stdmean.ps <- function(alpha, d, cor, w) {
 z <- qnorm(1 - alpha/2)
 n <- ceiling(4*(d^2*(1 + cor^2)/4 + 2*(1 - cor))*(z/w)^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 return(out)	
}


#  size.ci.ratio.mean.ps ===================================================== 
#' Sample size for a paired-samples mean ratio confidence interval 
#'
#'
#' @description
#' Computes the sample size required to estimate a ratio of population means 
#' with desired confidence interval precision in a paired-samples design.
#' Set the correlation planning value to the smallest value within a plausible
#' range for a conservatively large sample size. This function requires 
#' planning values for each mean and the sample size requirement is very 
#' sensitive to these planning values. Set the variance planning value to   
#' the largest value within a plausible range for a conservatively large  
#' sample size. 
#'
#'
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  var    planning value of average variance of the two measurements
#' @param  m1     planning value of mean for measurement 1
#' @param  m2     planning value of mean for measurement 2
#' @param  cor    planning value for correlation between measurements
#' @param  r      desired upper to lower confidence interval endpoint ratio
#'
#'
#' @return 
#' Returns the required sample size 
#'
#'
#' @examples
#' size.ci.ratio.mean.ps(.05, 400, 150, 100, .7, 1.2)
#'
#' # Should return:
#' #      Sample size
#' # [1,]          21
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.ci.ratio.mean.ps <- function(alpha, var, m1, m2, cor, r) {
 z <- qnorm(1 - alpha/2)
 n <- ceiling(8*var*(1/m1^2 + 1/m2^2 - 2*cor/(m1*m2))*(z/log(r))^2 + z^2/2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 return(out)
}


#  size.ci.lc.mean.ws ========================================================
#' Sample size for a within-subjects mean linear contrast confidence interval
#'
#'
#' @description
#' Computes the sample size required to estimate a linear contrast of 
#' population means with desired confidence interval precision in a 
#' within-subjects design. Set the variance planning value to the  
#' largest value within a plausible range for a conservatively large  
#' sample size. Set the correlation planning value to the smallest value
#' within a plausible range for a conservatively large sample size. 
#'
#'
#' @param  alpha  alpha level for 1-alpha confidence 
#' @param  var    planning value of average variance of the measurements  
#' @param  cor    planning value of average correlation between measurements
#' @param  w      desired confidence interval width
#' @param  q      vector of within-subjects contrast coefficients 
#'
#'
#' @return 
#' Returns the required sample size 
#'
#'
#' @examples
#' q <- c(.5, .5, -.5, -.5)
#' size.ci.lc.mean.ws(.05, 265, .8, 10, q)
#'
#' # Should return:
#' #      Sample size
#' # [1,]          11
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.ci.lc.mean.ws <- function(alpha, var, cor, w, q) {
 z <- qnorm(1 - alpha/2)
 k <- length(q)
 n <- ceiling(4*(1 - cor)*var*(t(q)%*%q)*(z/w)^2 + z^2/2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 return(out)
}


#  size.ci.lc.stdmean.ws =======================================================
#' Sample size for a within-subjects standardized mean linear contrast 
#' confidence interval
#'
#'
#' @description
#' Computes the sample size required to estimate a standardized linear 
#' contrast of population means with desired confidence interval precision
#' in a within-subjects design. Set the standardized mean difference planning 
#' value to the largest value within a plausible range for a conservatively 
#' large sample size. Set the correlation planning value to the smallest value
#' within a plausible range for a conservatively large sample size.
#'
#'
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  d      planning value of standardized linear contrast  
#' @param  cor    planning value of average correlation between measurements
#' @param  w      desired confidence interval width
#' @param  q      vector of within-subjects contrast coefficients 
#'
#'
#' @references
#' \insertRef{Bonett2009}{statpsych}
#'
#'
#' @return 
#' Returns the required sample size 
#'
#'
#' @examples
#' q <- c(.5, .5, -.5, -.5)
#' size.ci.lc.stdmean.ws(.05, 1, .7, .6, q)
#'
#' # Should return:
#' #      Sample size
#' # [1,]          26
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.ci.lc.stdmean.ws <- function(alpha, d, cor, w, q) {
 z <- qnorm(1 - alpha/2)
 a <- length(q)
 n <- ceiling(4*(d^2*(1 + (a - 1)*cor^2)/(2*a) + (1 - cor)*(t(q)%*%q))*(z/w)^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 return(out)
}


#  size.ci.cronbach ========================================================
#' Sample size for a Cronbach reliability confidence interval
#'
#'
#' Computes the sample size required to estimate a Cronbach reliability
#' with desired confidence interval precision.
#'
#'
#' @param  alpha  alpha value for 1-alpha confidence 
#' @param  rel    reliability planning value
#' @param  a      number of measurements
#' @param  w      desired confidence interval width
#'
#'
#' @references
#' \insertRef{Bonett2015}{statpsych}
#'
#'
#' @return 
#' Returns the required sample size
#'
#'
#' @examples
#' size.ci.cronbach(.05, .85, 5, .1)
#'
#' # Should return:
#' #      Sample size
#' # [1,]          89
#'  
#' 
#' @importFrom stats qf
#' @importFrom stats qnorm
#' @export
size.ci.cronbach <- function(alpha, rel, a, w) {
 z <- qnorm(1 - alpha/2)
 n0 <- ceiling((8*a/(a - 1))*(1 - rel)^2*(z/w)^2 + 2)
 df1 = n0 - 1
 df2 = n0*(a - 1)
 f1 = qf(1 - alpha/2, df1, df2)
 f2 = qf(1 - alpha/2, df2, df1)
 f0 <- 1/(1 - rel)
 ll <- 1 - f1/f0
 ul <- 1 - 1/(f0*f2)
 w0 <- ul - ll
 n <- ceiling((n0 - 2)*(w0/w)^2 + 2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 return(out)
}


#  size.ci.second =============================================================
#' Sample size for a second-stage confidence interval
#'
#'
#' @description
#' Computes the second-stage sample size required to obtain desired confidence 
#' interval precision. This function can use either the total sample size for 
#' all groups in the first stage sample or a single group sample size in the
#' first stage sample. If the total first-stage sample size is given, then the
#' function computes the total sample size required in the second-stage sample.
#' If a single group first-stage sample size is given, then the function 
#' computes the single-group sample size required in the second-stage sample. 
#' The second-stage sample is combined with the first-stage sample to obtain 
#' the desired confidence interval width.
#'
#'
#' @param  n0     first-stage sample size 
#' @param  w0     confidence interval width in first-stage sample
#' @param  w      desired confidence interval width
#'
#'
#' @return 
#' Returns the required sample size for the second-stage sample
#'
#'
#' @examples
#' size.ci.second(20, 5.3, 2.5)
#'
#' # Should return:
#' #      Second-stage sample size
#' # [1,]                       70
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.ci.second <- function(n0, w0, w) {
 n <- ceiling(((w0/w)^2 - 1)*n0)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 return(out)
}


# ==================== File 1: Sample Size for Desired Power ==================
#  size.test.mean1 ============================================================
#' Sample size for a test of a single mean
#'
#'
#' @description
#' Computes the sample size required to test a single population mean with 
#' desired power in a 1-group design. Set the variance planning value to the  
#' largest value within a plausible range for a conservatively large sample 
#' size.
#'
#'
#' @param  alpha  alpha level for hypothesis test 
#' @param  pow    desired power
#' @param  var    planning value of response variable variance
#' @param  es     planning value of mean minus null hypothesis value
#'
#'
#' @return 
#' Returns the required sample size 
#'
#'
#' @examples
#' size.test.mean1(.05, .9, 80.5, 7)
#'
#' # Should return:
#' #      Sample size
#' # [1,]          20
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.test.mean1 <- function(alpha, pow, var, es) {
 za <- qnorm(1 - alpha/2)
 zb <- qnorm(pow)
 n <- ceiling(var*(za + zb)^2/es^2 + za^2/2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 return(out)
}


#  size.test.mean2 ============================================================ 
#' Sample size for a test of a 2-group mean difference
#'
#'
#' @description
#' Computes the sample size in each group requires to test a difference in 
#' population means with desired power in a 2-group design. Set R =1 for equal 
#' sample sizes. Set the variance planning value to the largest value within a 
#' plausible range for a conservatively large sample size.
#'
#'
#' @param  alpha  alpha level for hypothesis test
#' @param  pow    desired power 
#' @param  var    planning value of average within-group variance
#' @param  es     planning value of mean difference
#' @param  R      n2/n1 ratio
#'
#'
#' @return 
#' Returns the required sample size per group
#'
#'	
#' @examples
#' size.test.mean2(.05, .95, 100, 10, 1) 
#'
#' # Should return:
#' #      n1  n2
#' # [1,] 27  27
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.test.mean2 <- function(alpha, pow, var, es, R) {
 za <- qnorm(1 - alpha/2)
 zb <- qnorm(pow)
 n1 <- ceiling(var*(1 + 1/R)*(za + zb)^2/es^2 + za^2/4)
 n2 <- ceiling(R*n1)
 out <- t(c(n1, n2))
 colnames(out) <- c("n1", "n2")
 return(out)
}


#  size.test.lc.mean.bs ====================================================== 
#' Sample size for a test of a between-subjects mean linear contrast
#'
#'
#' @description
#' Computes the sample size in each group (assuming equal sample sizes)
#' required to test a linear contrast of population means with desired 
#' power in a between-subjects design. Set the variance planning value
#' to the largest value within a plausible range for a conservatively 
#' large sample size.
#'
#'
#' @param  alpha  alpha level for hypothesis test 
#' @param  pow    desired power
#' @param  var    planning value of average within-group variance
#' @param  es     planning value of linear contrast of means
#' @param  v      vector of between-subjects contrast coefficients 
#'
#'
#' @return 
#' Returns the required sample size per group
#'
#'
#' @examples
#' v <- c(1, -1, -1, 1)
#' size.test.lc.mean.bs(.05, .90, 27.5, 5, v)
#'
#' # Should return:
#' #      Sample size per group
#' # [1,]                    47
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.test.lc.mean.bs <- function(alpha, pow, var, es, v) {
 za <- qnorm(1 - alpha/2)
 zb <- qnorm(pow)
 m <- length(v) - sum(v == 0)
 n <- ceiling(var*(t(v)%*%v)*(za + zb)^2/es^2 + za^2/(2*m))
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size per group"
 return(out)
}


#  size.equiv.mean2 ========================================================== 
#' Sample size for a 2-group mean equivalence test
#'
#'
#' @description
#' Computes the sample size in each group (assuming equal sample sizes) 
#' required to perform an equivalence test for the difference in population
#' means with desired power in a 2-group design. The value of h specifies a 
#' range of practical equivalence, -h to h, for the difference in population 
#' means. The planning value for the absolute mean difference must be less 
#' than h.  Equivalence tests often require a very large sample size. 
#' Equivalence tests usually use 2 x alpha rather than alpha (e.g., use 
#' alpha = .10 rather alpha = .05). Set the variance planning value to the
#' largest value within a plausible range for a conservatively large sample 
#' size.
#'
#'
#' @param   alpha  alpha level for hypothesis test 
#' @param   pow    desired power
#' @param   var    planning value of average within-group variance
#' @param   es     planning value of mean difference
#' @param   h      upper limit for range of practical equivalence
#'
#'
#' @return 
#' Returns the required sample size per group
#'
#'
#' @examples
#' size.equiv.mean2(.10, .80, 15, 2, 4)
#'
#' # Should return:
#' #      Sample size per group 
#' # [1,]                    50
#'  
#' 
#' @importFrom stats qnorm
#' @export
 size.equiv.mean2 <- function(alpha, pow, var, es, h) {
 za <- qnorm(1 - alpha)
 zb <- qnorm(1 - (1 - pow)/2)
 n <- ceiling(var*2*(za + zb)^2/(abs(h) - abs(es))^2 + za^2/4)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size per group"
 return(out)
}
 

#  size.supinf.mean2 ========================================================= 
#' Sample size for a 2-group mean superiority or noninferiority test
#'
#'
#' @description
#' Computes the sample size in each group (assuming equal sample sizes)  
#' required to perform a superiority or noninferiority test for the difference 
#' in population means with desired power in a 2-group design. For a 
#' superiority test, specify the upper limit (h) for the range of practical
#' equivalence and specify an effect size (es) such that es > h. For a 
#' noninferiority test, specify the lower limit (-h) for the range of 
#' practical equivalence and specify an effect size such that es > -h.  
#' Set the variance planning value to the largest value within a plausible 
#' range for a conservatively large sample size.
#'
#'
#' @param   alpha  alpha level for hypothesis test 
#' @param   pow    desired power
#' @param   var    planning value of average within-group variance
#' @param   es     planning value of mean difference
#' @param   h      upper or lower limit for range of practical equivalence
#'
#'
#' @return 
#' Returns the required sample size per group
#'
#'
#' @examples
#' size.supinf.mean2(.05, .80, 225, 9, 4)
#'
#' # Should return:
#' #      Sample size per group 
#' # [1,]                   143
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.supinf.mean2 <- function(alpha, pow, var, es, h) {
 za <- qnorm(1 - alpha/2)
 zb <- qnorm(pow)
 n <- ceiling(var*2*(za + zb)^2/(es - h)^2 + za^2/4)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 return(out)
}


#  size.test.mean.ps ========================================================== 
#' Sample size for a test of a paired-samples mean difference
#'
#'
#' @description
#' Computes the sample size required to test a difference in population means
#' with desired power in a paired-samples design. Set the correlation planning 
#' value to the smallest value within a plausible range for a conservatively 
#' large sample size. Set the variance planning value to the largest value 
#' within a plausible range for a conservatively large sample size.
#'
#'
#' @param  alpha  alpha level for hypothesis test 
#' @param  pow    desired power
#' @param  var    planning value of average variance of the two measurements
#' @param  es     planning value of mean difference 
#' @param  cor    planning value of correlation
#'
#'
#' @return 
#' Returns the required sample size 
#'
#'
#' @examples
#' size.test.mean.ps(.05, .80, 1.25, .5, .75) 
#'
#' # Should return:
#' #      Sample size
#' # [1,]          22
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.test.mean.ps <- function(alpha, pow, var, es, cor) {
 za <- qnorm(1 - alpha/2)
 zb <- qnorm(pow)
 n <- ceiling(2*var*(1 - cor)*(za + zb)^2/es^2 + za^2/2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 return(out)
}


#  size.test.lc.mean.ws ======================================================= 
#' Sample size for a test of a within-subjects mean linear contrast
#'
#'
#' @description
#' Computes the sample size required to test a linear contrast of population 
#' means with desired power in a within-subjects design. Set the average variance 
#' planning value to the largest value within a plausible range for a 
#' conservatively large sample size. Set the average correlation planning value to
#' the smallest value within a plausible range for a conservatively large sample 
#' size. 
#'
#'
#' @param  alpha  alpha level for hypothesis test 
#' @param  pow    desired power
#' @param  var    planning value of average variance of measurements  
#' @param  es     planning value of linear contrast of means
#' @param  cor    planning value of average correlation between measurements
#' @param  q      vector of with-subjects contrast coefficients 
#'
#'
#' @return 
#' Returns the required sample size 
#'
#'
#' @examples
#' q <- c(.5, .5, -.5, -.5)
#' size.test.lc.mean.ws(.05, .90, 50.7, 2, .8, q)
#'
#' # Should return:
#' #      Sample size
#' # [1,]          29
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.test.lc.mean.ws <- function(alpha, pow, var, es, cor, q) {
 za <- qnorm(1 - alpha/2)
 zb <- qnorm(pow)
 n <- ceiling(var*(1 - cor)*(t(q)%*%q)*(za + zb)^2/es^2 + za^2/2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 return(out)
}
	
	
#  size.equiv.mean.ps ========================================================== 
#' Sample size for a paired-samples mean equivalence test
#'
#'
#' @description
#' Computes the sample size required to perform an equivalence test for the
#' difference in population means with desired power in a paired-samples 
#' design. The value of h specifies a range of practical equivalence, -h to h, 
#' for the difference in population means. The planning value for the absolute
#' mean difference must be less than h. Equivalence tests often require a 
#' very large sample size. Equivalence tests usually use 2 x alpha rather than
#' alpha (e.g., use alpha = .10 rather alpha = .05). Set the correlation 
#' value to the smallest value within a plausible range for a conservatively 
#' large sample size. Set the variance planning value to the largest value 
#' within a plausible range for a conservatively large sample size.
#'
#'
#' @param   alpha  alpha level for hypothesis test 
#' @param   pow    desired power
#' @param   var    planning value of average within-group variance
#' @param   es     planning value of mean difference
#' @param   cor    planning value of the correlation between measurements
#' @param   h      upper limit for range of practical equivalence
#'
#'
#' @return 
#' Returns the required sample size
#'
#'
#' @examples
#' size.equiv.mean.ps(.10, .85, 15, .5, .7, 1.5)
#'
#' # Should return:
#' #      Sample size
#' # [1,]          68
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.equiv.mean.ps <- function(alpha, pow, var, es, cor, h) {
 za <- qnorm(1 - alpha)
 zb <- qnorm(1 - (1 - pow)/2)
 n <- ceiling(2*var*(1 - cor)*(za + zb)^2/(h - abs(es))^2 + za^2/2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 return(out)
}


#  size.supinf.mean.ps ======================================================== 
#' Sample size for a paired-samples mean superiority or noninferiority test
#'
#'
#' @description
#' Computes the sample size required to perform a superiority or noninferiority 
#' test for the difference in population means with desired power in a
#' paired-samples design. For a superiority test, specify the upper limit (h)
#' for the range of practical equivalence and specify an effect size (es) such
#' that es > h. For a noninferiority test, specify the lower limit (-h) for 
#' the range of practical equivalence and specify an effect size such that 
#' es > -h.  Set the correlation planning value to the smallest value within a
#' plausible range for a conservatively large sample size. Set the variance 
#' planning value to the largest value within a plausible range for a 
#' conservatively large sample size.
#'
#'
#' @param   alpha  alpha level for hypothesis test 
#' @param   pow    desired power
#' @param   var    planning value of average within-group variance
#' @param   es     planning value of mean difference
#' @param   cor    planning value of the correlation between measurements
#' @param   h      upper or lower limit for range of practical equivalence
#'
#'
#' @return 
#' Returns the required sample size
#'
#'
#' @examples
#' size.supinf.mean.ps(.05, .80, 225, 9, .75, 4)
#'
#' # Should return:
#' #      Sample size
#' # [1,]          38
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.supinf.mean.ps <- function(alpha, pow, var, es, cor, h) {
 za <- qnorm(1 - alpha/2)
 zb <- qnorm(pow)
 n <- ceiling(2*var*(1 - cor)*(za + zb)^2/(es - h)^2 + za^2/2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 return(out)
}


#  size.test.mann ============================================================ 
#' Sample size for a Mann-Whitney test
#'
#'
#' @description
#' Computes the sample size in each group (assuming equal sample sizes) 
#' required for the Mann-Whitney test with desired power. A planning value
#' of the Mann-Whitney parameter is required. In a 2-group experiment, this
#' parameter is the proportion of members in the population with scores that
#' would be larger under treatment 1 than treatment 2. In a 2-group 
#' nonexperiment where participants are sampled from two subpopulations of 
#' sizes N1 and N2, the parameter is the proportion of all N1 x N2 pairs in
#' which a member from subpopulation 1 has a larger score than a member from 
#' subpopulation 2.
#'
#'
#' @param  alpha  alpha level for hypothesis test 
#' @param  pow    desired power
#' @param  p      planning value of Mann-Whitney parameter
#'
#'
#' @references
#' \insertRef{Noether1987}{statpsych}
#'
#'
#' @return 
#' Returns the required sample size per group 
#'
#'
#' @examples
#' size.test.mann(.05, .90, .3)
#'
#' # Should return:
#' #      Sample size per group
#' # [1,]                    44
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.test.mann <- function(alpha, pow, p) {
 za <- qnorm(1 - alpha/2)
 zb <- qnorm(pow)
 es <- p - .5;
 n <- ceiling((za + zb)^2/(6*es^2))
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size per group"
 return(out)
}


#  size.test.sign1 =========================================================== 
#' Sample size for a 1-sample Sign test
#'
#'
#' @description
#' Computes the sample size required for a Sign test with desired power in a 
#' 1-sample design. A planning value of the 1-sample Sign test parameter value
#' is required. This parameter is the proportion of members in the population 
#' with scores greater than the hypothesized value. 
#'
#'
#' @param  alpha  alpha level for hypothesis test 
#' @param  pow    desired power
#' @param  p      planning value of proportion 
#'
#'
#' @return 
#' Returns the required sample size 
#'
#'
#' @examples
#' size.test.sign1(.05, .90, .3)
#'
#' # Should return:
#' #      Sample size
#' # [1,]          56
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.test.sign1 <- function(alpha, pow, p) {
 za <- qnorm(1 - alpha/2)
 zb <- qnorm(pow)
 es <- p - .5;
 n <- ceiling(p*(1 - p)*(za + zb)^2/es^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 return(out)
}


#  size.test.sign.ps ========================================================= 
#' Sample size for a paired-samples Sign test
#'
#'
#' @description
#' Computes sample size required for a Sign test with desired power in a 
#' paired-samples design. A planning value of the paired-samples Sign test 
#' parameter is required. In a paired-samples experiment, this parameter 
#' is the proportion of members in the population with scores that would be 
#' larger under treatment 1 than treatment 2. In a paired-samples 
#' nonexperiment, this parameter is the proportion of members in the  
#' population with measurement 1 scores that are larger than their
#' measurement 2 scores.
#'
#'
#' @param  alpha  alpha level for hypothesis test 
#' @param  pow    desired power
#' @param  p      planning value of proportion 
#'
#'
#' @return 
#' Returns the required sample size 
#'
#'
#' @examples
#' size.test.sign.ps(.05, .90, .75)
#'
#' # Should return:
#' #      Sample size
#' # [1,]          32
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.test.sign.ps <- function(alpha, pow, p) {
 za <- qnorm(1 - alpha/2)
 zb <- qnorm(pow)
 es <- p - .5;
 n <- ceiling(p*(1 - p)*(za + zb)^2/es^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 return(out)
}


#  size.test.cronbach ========================================================
#' Sample size to test a Cronbach reliability
#'
#'
#' Computes the sample size required to test a Cronbach reliability with
#' desired power. 
#'
#'
#' @param  alpha  alpha level for hypothesis test 
#' @param  pow    desired power
#' @param  rel    reliability planning value
#' @param  a      number of measurements
#' @param  h      null hypothesis value of reliability
#'
#'
#' @return 
#' Returns the required sample size
#'
#'
#' @references
#' \insertRef{Bonett2015}{statpsych}
#'
#'
#' @examples
#' size.test.cronbach(.05, .85, .80, 5, .7)
#'
#' # Should return:
#' #      Sample size
#' # [1,]         139
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.test.cronbach <- function(alpha, pow, rel, a, h) {
 za <- qnorm(1 - alpha/2)
 zb <- qnorm(pow)
 e <- (1 - rel)/(1 - h)
 n <- ceiling((2*a/(a - 1))*(za + zb)^2/log(e)^2 + 2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 return(out)
}


# ==========================  File 1: Miscellaneous ==========================
#  pi.score1 ================================================================= 
#' Prediction interval for one score
#'
#'
#' @description
#' Computes a prediction interval for the response variable score of one 
#' randomly selected member from the study population.
#'
#'
#' @param  alpha  alpha level for 1-alpha confidence 
#' @param  m      sample mean
#' @param  sd     sample standard deviation
#' @param  n      sample size
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Predicted - predicted score
#' * df - degrees of freedom
#' * LL - lower limit of the prediction interval
#' * UL - upper limit of the prediction interval
#'
#'
#' @examples
#' pi.score1(.05, 24.5, 3.65, 40)
#'
#' # Should return:
#' #      Predicted  df       LL       UL
#' # [1,]      24.5  39 17.02546 31.97454
#'  
#' 
#' @importFrom stats qt
#' @export
pi.score1 <- function(alpha, m, sd, n) {
 df <- n - 1
 tcrit <- qt(1 - alpha/2, df)
 se <- sqrt(sd^2 + sd^2/n)
 ll <- m - tcrit*se
 ul <- m + tcrit*se
 out <- t(c(m, df, ll, ul))
 colnames(out) <- c("Predicted", "df", "LL", "UL")
 return(out)
}


#  pi.score2 ================================================================== 
#' Prediction interval for a difference of scores
#'
#'
#' @description
#' For a 2-group experimental design, this function computes a prediction 
#' interval for how the response variable score for one randomly selected
#' person from the study population would differ under the two treatment
#' conditions. Both equal variance and unequal variance prediction intervals 
#' are computed.
#'
#'
#' @param  alpha  alpha level for 1-alpha confidence 
#' @param  m1     sample mean for group 1
#' @param  m2     sample mean for group 1
#' @param  sd1    sample standard deviation for group 1
#' @param  sd2    sample standard deviation for group 2
#' @param  n1     sample size for group 1
#' @param  n2     sample size for group 2
#'
#'
#' @return 
#' Returns a 2-row matrix. The columns are:
#' * Predicted - predicted difference in scores
#' * df - degrees of freedom
#' * LL - lower limit of the prediction interval
#' * UL - upper limit of the prediction interval
#'
#'
#' @references
#' \insertRef{Hahn1977}{statpsych}
#'
#'
#' @examples
#' pi.score2(.05, 29.57, 18.35, 2.68, 1.92, 40, 45)
#'
#' # Should return:
#' #                              Predicted       df       LL       UL
#' # Equal Variances Assumed:         11.22 83.00000 4.650454 17.78955
#' # Equal Variances Not Assumed:     11.22 72.34319 4.603642 17.83636
#'  
#' 
#' @importFrom stats qt
#' @export
pi.score2 <- function(alpha, m1, m2, sd1, sd2, n1, n2) {
 df1 <- n1 + n2 - 2
 est <- m1 - m2
 vp <- ((n1 - 1)*sd1^2 + (n2 - 1)*sd2^2)/df1
 se1 <- sqrt(2*vp + vp/n1 + vp/n2)
 tcrit1 <- qt(1 - alpha/2, df1)
 ll1 <- est - tcrit1*se1
 ul1 <- est + tcrit1*se1
 se2 <- sqrt(sd1^2 + sd2^2 + sd1^1/n1 + sd2^2/n2)
 c1 <- sd1^2 + sd1^2/n1
 c2 <- sd2^2 + sd2^2/n2
 df2 <- 1/((1/(n1 - 1))*(c1/(c1 + c2))^2 + (1/(n2 - 1))*(c2/(c1 + c2))^2)
 tcrit2 <- qt(1 - alpha/2, df2)
 ll2 <- est - tcrit2*se2
 ul2 <- est + tcrit2*se2
 out1 <- t(c(est, df1, ll1, ul1))
 out2 <- t(c(est, df2, ll2, ul2))
 out <- rbind(out1, out2)
 colnames(out) <- c("Predicted", "df", "LL", "UL")
 rownames(out) <- c("Equal Variances Assumed:", "Equal Variances Not Assumed:")
 return(out) 
}


#  random.sample ==============================================================
#' Generate a random sample
#'
#'
#' @description
#' Generates a random sample of participant IDs.
#'
#'
#' @param  popsize     study population size
#' @param  samsize     sample size
#'
#'
#' @return 
#' Returns a vector of randomly generated participant IDs
#'
#'
#' @examples
#' random.sample(3000, 25)
#'
#' # Should return:
#' #  [1]   37   94  134  186  212  408  485  697  722  781  998 1055 
#' # [13] 1182 1224 1273 1335 1452 1552 1783 1817 2149 2188 2437 2850 2936
#'  
#' 
#' @export
random.sample <- function(popsize, samsize) {
 out <- sort(sample(popsize, samsize, replace = FALSE, prob = NULL))
 return(out)
}


#  randomize ==================================================================
#' Randomize a sample into groups
#'
#'
#' @description
#' Randomly assigns a sample of participants into k groups.
#'
#'
#' @param  n   vector of sample sizes per group
#'
#'
#' @return 
#' Returns a vector of randomly generated group assignments.
#'
#'  
#' @examples
#' n <- c(10, 10, 5)
#' randomize(n)
#'
#' # Should return:
#' # [1] 2 3 2 1 1 2 3 3 2 1 2 1 3 1 1 2 3 1 1 2 2 1 1 2 2
#'  
#' 
#' @export
randomize <- function(n) {
 k <- length(n)
 out <- sample(rep(1:k, n))
 return(out)
}


#  random.y =================================================================== 
#' Generate random sample of scores
#'
#'
#' @description
#' Generates a random sample of scores from a normal distribution with a
#' specified population mean and standard deviation.
#'
#'
#' @param  n     sample size
#' @param  m     population mean of scores
#' @param  sd    population standard deviation of scores
#' @param  min   minimum allowable value
#' @param  max   maximum allowable value
#' @param  dec   number of decimal points 
#'
#'
#' @return 
#' Returns a vector of randomly generated scores.
#'
#'  
#' @examples
#' random.y(10, 3.6, 2.8, 1, 7, 0) 
#'
#' # Should return:
#' # [1] 2 7 7 1 6 3 1 3 2 1
#'  
#' 
#' @importFrom stats rnorm
#' @export
random.y <- function(n, m, sd, min, max, dec) {
 y <- sd*rnorm(n, 0, 1) + m
 y1 <- 1 - as.integer(y < min)
 y <-  y*y1 + (1 - y1)*1
 y2 <- 1 - as.integer(y > max)
 y = y*y2 + (1 - y2)*max
 out <- round(y, dec)
 return(out)
} 


#  ci.var.upper =============================================================== 
#' Upper confidence limit of variance
#'
#'
#' @description
#' Computes an upper limit for a population variance using a sample variance 
#' from a sample of size n in a prior study. The upper limit can be used as 
#' a variance planning value in sample size functions that require a variance 
#' planning value.
#'
#'
#' @param  alpha  alpha value for 1-alpha confidence (one-sided)
#' @param  var    sample variance
#' @param  n      sample size

#' @return 
#' Returns an upper limit (UL) variance planning value
#'
#'
#' @examples
#' ci.var.upper(.25, 15, 60)
#'
#' # Should return:
#' #            UL
#' # [1,] 17.23264
#'  
#' 
#' @importFrom stats qchisq
#' @export
ci.var.upper <- function(alpha, var, n) {
 chi <- qchisq(alpha, (n - 1))
 ul <- (n - 1)*var/chi
 out <- matrix(ul, nrow = 1, ncol = 1)
 colnames(out) <- "UL"
 return(out)
}


#  etasqr.adj =================================================================
#' Bias adjusts an eta-squared estimate
#'
#'
#' @description
#' Computes an approximate bias adjustment for eta-squared. This adjustment
#' can be applied to eta-squared, partial-eta squared, and generalized
#' eta-squared estimates.
#'
#'
#' @param   etasqr    unadjusted eta-square estimate
#' @param   dfeffect  degrees of freedom for the effect
#' @param   dferror   error degrees of freedom
#'
#'
#' @return 
#' Returns a bias adjusted eta-squared estimate
#'
#'
#' @examples
#' etasqr.adj(.315, 2, 42)
#'
#' # Should return:
#' #      Adjusted eta-squared
#' # [1,]             0.282381
#'  
#' 
#' @export
etasqr.adj <- function(etasqr, dfeffect, dferror) {
 adj <- 1 - (dferror + dfeffect)*(1 - etasqr)/dferror
 out <- matrix(adj, nrow = 1, ncol = 1)
 colnames(out) <- "Adjusted eta-squared"
 return(out)
}


#  test.anova1.bs =============================================================
#' Between-subjects F statistic and eta-squared from summary information 
#'
#'
#' @description
#' Computes the F statistic, p-value, eta-squared, and adjusted eta-squared 
#' for the main effect of Factor A in a one-way between-subjects ANOVA using
#' the sample means, sample standard deviations, and sample sizes.  
#'
#'
#' @param   m       vector of sample means
#' @param   sd      vector of sample standard deviations
#' @param   n       vector of sample sizes
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * F - F statistic for test of null hypothesis
#' * dfA - degrees of freedom for between-subjects factor
#' * dfE - error degrees of freedom
#' * dfA - degrees of freedom for between-subjects factor
#' * p - p-value for F-test
#' * eta-squared - estimate of eta-squared
#' * adj eta-squared - a bias adjusted estimate of eta-squared
#'
#'
#' @examples
#' m <- c(12.4, 8.6, 10.5)
#' sd <- c(3.84, 3.12, 3.48)
#' n <- c(20, 20, 20)
#' test.anova1.bs(m, sd, n)
#'
#' #  Should return:
#' #             F dfA  dfE           p eta-squared adj eta-squared
#' # [1,] 5.919585   2   57 0.004614428   0.1719831       0.1429298
#'  
#' 
#' @importFrom stats pf
#' @export
test.anova1.bs <- function(m, sd, n) {
  a <- length(m)
  nt <- sum(n)
  dfe <- nt - a
  dfa <- a - 1
  v <- sd^2
  grandmean <- sum(n*m)/nt
  SSe <- sum((n - 1)*v)
  MSe <- SSe/dfe
  SSa <- sum(n*(m - grandmean)^2)
  MSa <- SSa/dfa
  F <- MSa/MSe
  p <- 1 - pf(F, dfa, dfe)
  etasqr <- SSa/(SSa + SSe)
  adjetasqr <- 1 - (dfa + dfe)*(1 - etasqr)/dfe
  out <- t(c(F, dfa, dfe, p, etasqr, adjetasqr))
  colnames(out) <- c("F", "dfA",  "dfE", "p", "eta-squared", "adj eta-squared")
  return(out)
}



 
 

 
