# =========================== Confidence Intervals ===========================
#  ci.mean1  =================================================================
#' Confidence interval for a single mean
#'
#'
#' @description
#' Computes a confidence interval for a population mean using the estimated 
#' mean, estimated standard deviation, and sample size. Use the t.test function
#' for raw data input.
#'
#'  
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  m	  estimated mean 
#' @param  sd	  estimated standard deviation
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
#' @param  m	  estimated mean 
#' @param  sd	  estimated standard deviation
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
#' Computes equal variance and unequal variance confidence intervals for a 
#' population 2-group mean difference using the estimated means, estimated 
#' standard deviations, and sample sizes. Use the t.test function for raw data
#' input.
#'
#'  
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  m1     estimated mean for group 1
#' @param  m2     estimated mean for group 2
#' @param  sd1    estimated standard deviation for group 1
#' @param  sd2    estimated standard deviation for group 2
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
#' of means. This function computes both unequal variance and equal
#' variance confidence intervals and test statistics. A Satterthwaite 
#' adjustment to the degrees of freedom is used with the unequal variance 
#' method. 
#'
#'
#' @param     alpha  	alpha level for 1-alpha confidence
#' @param     m     	vector of estimated group means
#' @param     sd    	vector of estimated group standard deviations
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
#' using estimated means, estimated standard deviations, and samples sizes as 
#' input. A Satterthwaite adjustment to the degrees of freedom is used to 
#' improve the accuracy of the confidence intervals. 
#'
#'
#' @param  alpha   alpha level for simultaneous 1-alpha confidence
#' @param  m       vector of estimated group means
#' @param  sd      vector of estimated group standard deviations
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
#' m <- c(12.86, 17.57, 26.29, 30.21)
#' sd <- c(13.185, 12.995, 14.773, 15.145)
#' n <- c(20, 20, 20, 20)
#' ci.tukey(.05, m, sd, n)
#'
#' # Should return:
#' #     Estimate       SE          t       df           p        LL         UL
#' # 1 2    -4.71 4.139530 -1.1378102 37.99200 0.668806358 -15.83085  6.4108517
#' # 1 3   -13.43 4.427673 -3.0331960 37.51894 0.021765570 -25.33172 -1.5282764
#' # 1 4   -17.35 4.490074 -3.8640790 37.29278 0.002333937 -29.42281 -5.2771918
#' # 2 3    -8.72 4.399497 -1.9820446 37.39179 0.212906199 -20.54783  3.1078269
#' # 2 4   -12.64 4.462292 -2.8326248 37.14275 0.035716267 -24.64034 -0.6396589
#' # 3 4    -3.92 4.730817 -0.8286096 37.97652 0.840551420 -16.62958  8.7895768
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
 Estimate <- (-1)*mean[lower.tri(mean)]
 v1 <- outer(v1, v1, "+")
 v2 <- outer(v2, v2, "+")
 df = v1^2/v2
 df <- df[lower.tri(df)]
 SE <- sqrt(v1[lower.tri(v1)])
 t <- Estimate/SE
 q <- qtukey(p = 1 - alpha, nmeans = a, df = df)/sqrt(2)
 p <- 1 - ptukey(sqrt(2)*abs(t), nmeans = a, df = df)
 LL <- Estimate - q*SE
 UL <- Estimate + q*SE
 pair <- t(combn(seq(1:a), 2))
 out <- cbind(pair, Estimate, SE, t, df, p, LL, UL)
 rownames(out) <- rep("", a*(a - 1)/2)
 return(out)
}


#  ci.slope.mean.bs ===========================================================
#' Confidence interval for the slope of means in a single-factor design with a
#' quantitative between-subjects factor
#' 
#' 
#' @description
#' Computes a test statistic and confidence interval for the slope of means in 
#' a single-factor design with a quantitative between-subjects factor. This 
#' function computes both the unequal variance and equal variance confidence
#' intervals and test statistics. A Satterthwaite adjustment to the degrees of
#' freedom is used with the unequal variance method. 
#'
#'
#' @param     alpha  	alpha level for 1-alpha confidence
#' @param     m     	vector of sample means
#' @param     sd    	vector of sample standard deviations
#' @param     n     	vector of sample sizes
#' @param     x     	vector of numeric predictor variable values
#' 
#'
#' @return 
#' Returns a 2-row matrix. The columns are:
#' * Estimate - estimated slope
#' * SE - standard error
#' * t - t test statistic
#' * df - degrees of freedom
#' * p - p-value
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' m <- c(33.5, 37.9, 38.0, 44.1)
#' sd <- c(3.84, 3.84, 3.65, 4.98)
#' n <- c(10,10,10,10)
#' x <- c(5, 10, 20, 30)
#' ci.slope.mean.bs(.05, m, sd, n, x)
#'
#' # Should return:
#' #                               Estimate         SE        t       df
#' # Equal Variances Assumed:     0.3664407 0.06770529 5.412290 36.00000
#' # Equal Variances Not Assumed: 0.3664407 0.07336289 4.994905 18.65826
#' #                                         p        LL        UL
#' # Equal Variances Assumed:     4.242080e-06 0.2291280 0.5037534
#' # Equal Variances Not Assumed: 8.468223e-05 0.2126998 0.5201815
#'
#'
#' @importFrom stats qt
#' @importFrom stats pt
#' @export
ci.slope.mean.bs <- function(alpha, m, sd, n, x) {
 mx <- mean(x)
 ssx <- sum((x - mx)^2)
 v <- (x - mx)/ssx
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
#' @param  m1     estimated mean for group 1
#' @param  m2     estimated mean for group 2
#' @param  sd1    estimated standard deviation for group 1
#' @param  sd2    estimated standard deviation for group 2
#' @param  n1     sample size for group 1
#' @param  n2     sample size for group 2
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
#' in a 2-group nonexperimental design with stratified random sampling (a random
#' sample of a specificied size from each subpopulation) using a square root 
#' weighted variance standardizer or single group standard deviation 
#' standardizer.  Equality of variances is not assumed.
#'
#'
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  m1     estimated mean for group 1
#' @param  m2     estimated mean for group 2
#' @param  sd1    estimated standard deviation for group 1
#' @param  sd2    estimated standard deviation for group 2
#' @param  n1     sample size for group 1
#' @param  n2     sample size for group 2
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
#' of means in a between-subjects design. The unweighted standardizer is 
#' recommended in experimental designs. The weighted standardizer is
#' recommended in nonexperimental designs with simple random sampling. The  
#' group 1 standardizer is useful in both experimental and nonexperimental
#' designs. Equality of variances is not assumed.
#'
#'
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  m      vector of estimated group means
#' @param  sd     vector of estimated group standard deviation
#' @param  n      vector of sample sizes
#' @param  v      vector of between-subjects contrast coefficients
#'
#'
#' @return 
#' Returns a 3-row matrix. The columns are:
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
#' # Group 1 standardizer:    -1.273810 0.4849842 -2.343781 -0.4426775
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
 n1 <- matrix(n, 1, a)[1,1]
 s1 <- matrix(sd, 1, a)[1,1]
 adj1 <- 1 - 3/(4*df - 1)
 adj2 <- 1 - 3/(4*n1 - 5)
 sp <- sqrt(sum((n - 1)*var)/df)
 est1 <- (t(v)%*%m)/s
 est2 <- (t(v)%*%m)/sp
 est3 <- (t(v)%*%m)/s1 
 est1u <- adj1*est1
 est2u <- adj1*est2
 est3u <- adj2*est3
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
 a1 <- est3^2/(2*n1 - 2)
 a2 <- sum((v^2*var/(n - 1)))/s1^2
 se3 <- sqrt(a1 + a2)
 ll3 <- est3 - z*se3
 ul3 <- est3 + z*se3
 out1 <- t(c(est1u, se1, ll1, ul1))
 out2 <- t(c(est2u, se2, ll2, ul2))
 out3 <- t(c(est3u, se3, ll3, ul3))
 out <- rbind(out1, out2, out3)
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 rownames(out) <- c("Unweighted standardizer:", "Weighted standardizer:", "Group 1 standardizer:")
 return(out)
}


#  ci.mean.ps ================================================================
#' Confidence interval for a paired-samples mean difference
#'
#' @description
#' Computes a confidence interval for a population paired-samples mean 
#' difference using the estimated means, estimated standard deviations, 
#' estimated correlation, and sample size. Use the t.test function for 
#' raw data input.
#'
#'  
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  m1     estimated mean for measurement 1
#' @param  m2     estimated mean for measurement 2
#' @param  sd1    estimated standard deviation for measurement 1
#' @param  sd2    estimated standard deviation for measurement 2
#' @param  cor    estimated correlation between measurements
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
#' @param  m1     estimated mean of measurement 1
#' @param  m2     estimated mean of measurement 2
#' @param  sd1    estimated standard deviation of measurement 1
#' @param  sd2    estimated standard deviation of measurement 2
#' @param  cor    estimated correlation between measurements
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
#' Computes confidence intervals for two types of population standardized 
#' linear contrast of means (unweighted standardizer and level 1 standardizer)
#' in a within-subjects design. Equality of variances is not assumed, but the
#' correlations among the repeated measures are assumed to be approximately 
#' equal.
#'
#'
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  m      vector of within-subjects estimated means
#' @param  sd     vector of within-subjects estimated standard deviations
#' @param  cor    average estimated correlation of all measurement pairs
#' @param  n      sample size
#' @param  q      vector of within-subjects contrast coefficients
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
#' q <- c(.5, .5, -.5, -.5)
#' ci.lc.stdmean.ws(.05, m, sd, .672, 20, q)
#'
#' # Should return:
#' #                           Estimate        SE        LL         UL
#' # Unweighted standardizer: -1.266557 0.2096351 -1.712140 -0.8903860
#' # Level 1 standardizer:    -1.337500 0.2662156 -1.915002 -0.8714561
#'
#'
#' @importFrom stats qnorm
#' @export
ci.lc.stdmean.ws <- function(alpha, m, sd, cor, n, q) {
 z <- qnorm(1 - alpha/2)
 a <- length(m)
 df <- n - 1
 s1 <- matrix(sd, 1, a)[1,1]
 adj1 <- sqrt((n - 2)/df)
 adj2 <- 1 - 3/(4*n - 5)
 s <- sqrt(sum(sd^2)/a)
 est1 <- (t(q)%*%m)/s
 est1u <- adj1*est1
 v1 <- est1^2/(2*a^2*s^4*df)
 v2 <- sum(sd^4)
 v0 <- sd^2%*%t(sd^2)
 v3 <- cor^2*sum((v0 - diag(diag(v0))))
 # error in v3 formula corrected in version 1.3
 v4 <- sum(q^2*sd^2)
 v5 <- cor*t(q*sd)%*%(q*sd)
 se1 <- sqrt(v1*(v2 + v3) + (v4 - v5)/(df*s^2))
 ll1 <- est1 - z*se1
 ul1 <- est1 + z*se1
 est2 <- (t(q)%*%m)/s1
 est2u <- adj2*est2
 v1 <- est2^2/(2*df)
 v4 <- sum(q^2*sd^2)
 v5 <- cor*t(q*sd)%*%(q*sd)
 se2 <- sqrt(v1 + (v4 - v5)/(df*s1^2))
 ll2 <- est2 - z*se2
 ul2 <- est2 + z*se2
 out1 <- t(c(est1u, se1, ll1, ul1))
 out2 <- t(c(est2u, se2, ll2, ul2))
 out <- rbind(out1, out2)
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 rownames(out) <- c("Unweighted standardizer:", "Level 1 standardizer:")
 return(out)
}


#  ci.mad1 ====================================================================
#' Confidence interval for a single mean absolute deviation 
#'
#'
#' @description
#' Computes a confidence interval for a population mean absolute deviation 
#' from the median (MAD). The MAD is a robust alternative to the standard 
#' deviation.
#'
#'
#' @param  alpha   alpha level for 1-alpha confidence
#' @param  y       vector of scores
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimated mean absolute deviation
#' * SE - standard error
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
#' #       Estimate       SE       LL       UL
#' # [1,]      12.5 2.876103 7.962667 19.62282
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
 se.mad <- c*mad*se
 out <- t(c(c*mad, se.mad, ll, ul))
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 return(out)
}


#  ci.ratio.mad2 ==============================================================
#' Confidence interval for a 2-group ratio of mean absolute deviations
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


#  ci.ratio.sd2 ================================================================
#' Confidence interval for a 2-group ratio of standard deviations 
#'
#'
#' @description
#' Computes a robust confidence interval for a ratio of population standard 
#' deviations in a 2-group design. This function is a modification of the
#' confidence interval proposed by Bonett (2006). The original Bonett method 
#' used a pooled kurtosis estimate in the standard error that assumed equal 
#' variances, which limited the confidence interval's use to tests of equal 
#' population variances and equivalence tests. This function uses a pooled 
#' kurtosis estimate that does not assume equal variances and provides a 
#' useful confidence interval for a ratio of standard deviations under general 
#' conditions. This function requires of minimum sample size of four per
#' group but sample sizes of at least 10 per group are recommended.
#'
#'
#' @param  alpha   alpha level for 1-alpha confidence
#' @param  y1      vector of scores for group 1
#' @param  y2      vector of scores for group 2
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * SD1 - estimated SD from group 1
#' * SD2 - estimated SD from group 2
#' * SD1/SD2 - estimate of SD ratio
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#'
#' @references
#' \insertRef{Bonett2006b}{statpsych}                   
#'
#'
#' @examples
#' y1 <- c(32, 39, 26, 35, 43, 27, 40, 37, 34, 29)
#' y2 <- c(36, 44, 47, 42, 49, 39, 46, 31, 33, 48)
#' ci.ratio.sd2(.05, y1, y2)
#'
#' # Should return:
#' #           SD1      SD2    SD1/SD2       LL       UL
#' # [1,] 5.711587 6.450667  0.8854257 0.486279 1.728396
#'
#'
#' @importFrom stats qnorm
#' @importFrom stats sd
#' @export
ci.ratio.sd2 <- function(alpha, y1, y2) {
  z <- qnorm(1 - alpha/2)
  sd1 <- sd(y1)
  sd2 <- sd(y2)
  v1 <- sd1^2
  v2 <- sd2^2
  n1 <- length(y1)
  n2 <- length(y2)
  if (min(n1, n2) < 5) {stop("sample size too small")}
  n <- n1 + n2 
  t1 <- 1/(2*sqrt(n1 - 4))
  t2 <- 1/(2*sqrt(n2 - 4))
  m1 <- mean(y1, trim = t1)
  m2 <- mean(y2, trim = t2)
  c0 <- (n1/(n1 - z))*((n2 - z)/n2)
  c1 <- (n1 - 3)/n1
  c2 <- (n2 - 3)/n2
  a1 <- sum((y1 - m1)^4)
  a2 <- sum((y2 - m2)^4)
  rL <- 1
  rU <- 1
  for(i in 1:10) {
    kurL <- n*(a1 + rL^4*a2)/((n1 - 1)*v1 + rL^2*(n2 - 1)*v2)^2
    kurU <- n*(a1 + rU^4*a2)/((n1 - 1)*v1 + rU^2*(n2 - 1)*v2)^2
    seL <- sqrt((kurL - c1)/(n1 - 1) + (kurL - c2)/(n2 - 1))
    seU <- sqrt((kurU - c1)/(n1 - 1) + (kurU - c2)/(n2 - 1))
    lr <- log(c0*v1/v2)
    ll <- sqrt(exp(lr - z*seL))
    ul <- sqrt(exp(lr + z*seU))
    rL <- ll
    rU <- ul
  }
 out <- t(c(sd1, sd2, sd1/sd2, ll, ul))
 colnames(out) <- c("SD1", "SD2", "SD1/SD2", "LL", "UL")
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
#' which is defined as a mean absolute deviation from the median divided by a
#' median. The coefficient of dispersion assumes ratio-scale scores and is a 
#' robust alternative to the coefficient of variation. An approximate
#' standard error is recovered from the confidence interval.
#'
#' @param  alpha   alpha level for 1-alpha confidence
#' @param  y       vector of scores
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimated coefficient of dispersion
#' * SE - recovered standard error
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
#' #        Estimate        SE        LL       UL
#' # [1,]  0.5921053 0.1814708 0.3813259 1.092679
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
 ll = exp(L1 - U2)
 ul = exp(U1 - L2)
 se <- (ul - ll)/(2*z)
 out = t(c(cod, se, ll, ul))
 colnames(out) = c("Estimate", "SE", "LL", "UL")
 return(out)
}


# ci.cod2 ====================================================================
#' Confidence interval for a ratio of dispersion coefficients in a 2-group
#' design
#'
#'           
#' @description
#' Computes a confidence interval for a ratio of population dispersion
#' coefficients (MAD/median) in a 2-group design. Ratio-scale scores are
#' assumed.
#'
#' @param  alpha   alpha level for 1-alpha confidence
#' @param  y1      vector of scores in gorup 1
#' @param  y2      vector of scores in gorup 2
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * COD1 - estimated coefficient of dispersion in group 1
#' * COD2 - estimated coefficient of dispersion in group 2
#' * COD1/COD2 - estimated ratio of dispersion coefficients
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @examples
#' y1 <- c(32, 39, 26, 35, 43, 27, 40, 37, 34, 29)
#' y2 <- c(36, 44, 47, 42, 49, 39, 46, 31, 33, 48)
#' ci.cod2(.05, y1, y2)
#'
#' # Should return:
#' #           COD1      COD2 COD1/COD2       LL       UL
#' # [1,] 0.1333333 0.1232558  1.081761 0.494964 2.282254
#'
#'
#' @importFrom stats qnorm
#' @export
ci.cod2 <-function(alpha, y1, y2) {
 z <- qnorm(1 - alpha/2)
 n1 <- length(y1)
 c1 <- n1/(n1 - 1)
 a1 <- round((n1 + 1)/2 - sqrt(n1))
 b1 <- n1 - a1 + 1
 med1 <- median(y1)
 m1 <- mean(y1)
 v1 <- var(y1)
 tau1 <- mean(abs(y1 - med1))
 del1 <- (m1 - med1)/tau1
 gam1 <- v1/tau1^2
 cod1 <- tau1/med1
 y1 <- sort(y1)
 vareta1 <- ((log(y1[a1]) - log(y1[b1]))/4)^2
 se11 <- sqrt(vareta1)
 vartau1 <- (gam1 + del1^2 - 1)/n1
 se21 <- sqrt(vartau1)
 cov1 <- del1*se11/sqrt(n1)
 k1 <- sqrt(vareta1 + vartau1 - 2*cov1)/(se11 + se21)
 a21 <- round((n1 + 1)/2 - k1*z*sqrt(n1/4))
 b21 <- n1 - a21 + 1
 L21 <- log(y1[a21])
 U21 <- log(y1[b21])
 L11 <- log(c1*tau1) - k1*z*se21
 U11 <- log(c1*tau1) + k1*z*se21
 LL1 <- L11 - U21
 UL1 <- U11 - L21
 n2 <- length(y2)
 c2 <- n2/(n2 - 1)
 a2 <- round((n2 + 1)/2 - sqrt(n2))
 b2 <- n2 - a2 + 1
 med2 <- median(y2)
 m2 <- mean(y2)
 v2 <- var(y2)
 tau2 <- mean(abs(y2 - med2))
 del2 <- (m2 - med2)/tau2
 gam2 <- v2/tau2^2
 cod2 <- tau2/med2
 y2 <- sort(y2)
 vareta2 <- ((log(y2[a1]) - log(y2[b1]))/4)^2
 se12 <- sqrt(vareta2)
 vartau2 <- (gam2 + del2^2 - 1)/n2
 se22 <- sqrt(vartau2)
 cov2 <- del2*se12/sqrt(n2)
 k2 <- sqrt(vareta2 + vartau2 - 2*cov2)/(se12 + se22)
 a22 <- round((n2 + 1)/2 - k2*z*sqrt(n2/4))
 b22 <- n2 - a22 + 1
 L22 <- log(y2[a22])
 U22 <- log(y2[b22])
 L12 <- log(c2*tau2) - k2*z*se22
 U12 <- log(c2*tau2) + k2*z*se22
 LL2 <- L12 - U22
 UL2 <- U12 - L22
 LL <- exp(log(cod1/cod2) - sqrt((log(cod1) - LL1)^2 + (log(cod2) - UL2)^2))
 UL <- exp(log(cod1/cod2) + sqrt((log(cod1) - UL1)^2 + (log(cod2) - LL2)^2))
 out <- t(c(cod1, cod2, cod1/cod2, LL, UL))
 colnames(out) <- c("COD1", "COD2", "COD1/COD2", "LL", "UL")
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
#' * Estimate - estimated median
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
#' #      Estimate       SE LL UL
#' # [1,]       20 4.270922 10 30
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
 a <- round(n/2 - sqrt(n))
 if (a < 1) {a = 1}
 ll1 <- y[a]
 ul1 <- y[n - a + 1]
 p <- pbinom(a - 1, size = n, prob = .5)
 z0 <- qnorm(1 - p)
 se <- (ul1 - ll1)/(2*z0)
 out <- t(c(median, se, ll, ul))
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 return(out)
}


#  ci.median2 =================================================================
#' Confidence interval for a 2-group median difference
#'
#'
#' @description
#' Computes a confidence interval for a difference of population medians in
#' a 2-group design.
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
#' between-subjects design using estimated medians and their standard errors.
#' The sample median and standard error for each group can be computed using
#' the ci.median1 function.
#'
#'
#' @param  alpha   alpha level for 1-alpha confidence
#' @param  m       vector of estimated group medians
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
#' for each median and the covariance between the two estimated medians.
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
#' * COV - covariance of the two `estimated medians
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
#' #        SE1      SE2      COV
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
 colnames(out) <- c("Median1", "Median2", "Median1-Median2", "SE", "LL", "UL", "SE1", "SE2", "COV")
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


#  ci.sign1 ================================================================== 
#' Confidence interval for the parameter of the one-sample sign test
#'
#'
#' @description
#' Computes adjusted Wald interval for the population proportion of 
#' quantitative scores that are greater than the null hypothesis value
#' of the population median in a one-sample sign test.
#'
#'
#' @param   alpha   alpha level for 1-alpha confidence
#' @param   y       vector of y scores
#' @param   h       null hypothesis value for population median
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - adjusted estimate of proportion
#' * SE - adjusted standard error
#' * LL - lower limit of adjusted Wald confidence interval
#' * UL - upper limit of adjusted Wald confidence interval
#'
#'
#' @references
#' Agresti, A, & Coull, BA (1998) Approximate is better than “exact” for 
#' interval estimation of binomial proportions. American Statistician, 52,
#' 119–126. doi: 10.1080/00031305.1998.10480550
#'
#'
#' @examples
#' y <- c(30, 20, 15, 10, 10, 60, 20, 25, 20, 30, 10, 5, 50, 40, 20, 10,
#'         0, 20, 50)
#' ci.sign1(.05, y, 9)
#'
#' # Should return:
#' #      Estimate        SE        LL        UL
#' # [1,] 0.826087 0.0790342 0.6711828 0.9809911
#'
#'
#' @importFrom stats qnorm
#' @export
ci.sign1 <- function(alpha, y, h) {
 z <- qnorm(1 - alpha/2)
 f <- sum(as.integer(y > h))
 n <- length(y)
 p.mle <- f/n
 se.mle <- sqrt(p.mle*(1 - p.mle)/n)
 p.adj <- (f + 2)/(n + 4)
 se.adj <- sqrt(p.adj*(1 - p.adj)/(n + 4))
 LL.adj <- p.adj - z*se.adj
 UL.adj <- p.adj + z*se.adj
 if (LL.adj < 0) {LL.adj = 0}
 if (UL.adj > 1) {UL.adj = 1}
 out <- t(c(p.adj, se.adj, LL.adj, UL.adj))
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 return(out)
}


#  ci.mann ====================================================================
#' Confidence interval for a Mann-Whitney parameter
#'
#'
#' @description
#' Computes a distribution-free confidence interval for the Mann-Whitney 
#' parameter (a "common language effect size"). In a 2-group experiment, this
#' parameter is the proportion of members in the population with scores that 
#' would be higher under treatment 1 than treatment 2. In a 2-group 
#' nonexperiment where participants are sampled from two subpopulations of 
#' sizes N1 and N2, the parameter is the proportion of all N1 x N2 pairs in 
#' which a member from subpopulation 1 has a larger score than a member from 
#' subpopulation 2.
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
#' y2 <- c(36, 44, 47, 42, 49, 39, 46, 31, 33, 48)
#' y1 <- c(32, 39, 26, 35, 43, 27, 40, 37, 34, 29)
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
 u <- sum(r1) - n1*(n1 + 1)/2
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
#' @param   m      vector of estimated group means 
#' @param   sd     vector of estimated group standard deviations 
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
#' @param  r      number of measurements (items, raters, etc.)
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
#' #               SE        LL        UL
#' # [1,]  0.02456518 0.7971254 0.8931436   
#'  
#' 
#' @importFrom stats qf
#' @export
ci.cronbach <- function(alpha, rel, r, n) {
 se <- sqrt((2*r*(1 - rel)^2)/((r - 1)*(n - 2)))
 df1 = n - 1
 df2 = n*(r - 1)
 f1 = qf(1 - alpha/2, df1, df2)
 f2 = qf(1 - alpha/2, df2, df1)
 f0 <- 1/(1 - rel)
 ll <- 1 - f1/f0
 ul <- 1 - 1/(f0*f2)
 out <- t(c(se, ll, ul))
 colnames(out) = c("SE", "LL", "UL")
 return(out)
}


#  ci.reliablity =============================================================
#' Confidence interval for a reliability coefficient
#'
#'
#' @description
#' Computes a confidence interval for a population reliability coefficient
#' such as Cronbach's alpha or McDonald's omega using an estimate of the
#' reliablity and its standard error. The standard error can be a robust
#' standard error or bootstrap standard error obtained from an SEM program.
#'
#'  
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  rel    estimated reliability  
#' @param  se     standard error of reliability
#' @param  n	  sample size
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' ci.reliability(.05, .88, .0147, 100)
#'
#' # Should return:
#' #             LL        UL
#' # [1,] 0.8489612 0.9065575
#'  
#' 
#' @importFrom stats qf
#' @export
ci.reliability <- function(alpha, rel, se, n) {
 z <- qnorm(1 - alpha/2)
 b <- log(n/(n - 1))
 ll <- 1 - exp(log(1 - rel) - b + z*sqrt(se^2/(1 - rel)^2))
 ul <- 1 - exp(log(1 - rel) - b - z*sqrt(se^2/(1 - rel)^2))
 out <- t(c(ll, ul))
 colnames(out) = c("LL", "UL")
 return(out)
}


#  ci.etasqr ================================================================= 
#' Confidence interval for eta-squared
#'
#'
#' @description
#' Computes a confidence interval for a population eta-squared, partial 
#' eta-squared, or generalized eta-squared in a fixed-factor between-subjects 
#' design. An approximate bias adjusted estimate is computed, and an
#' approximate standard error is recovered from the confidence interval.
#'
#'
#' @param  alpha    alpha value for 1-alpha confidence
#' @param  etasqr   estimated eta-squared
#' @param  df1      degrees of freedom for effect
#' @param  df2      error degrees of freedom
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Eta-squared - estimate of eta-squared 
#' * adj Eta-squared - bias adjusted eta-squared estimate
#' * SE - recovered standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @examples
#' ci.etasqr(.05, .241, 3, 116)
#'
#' # Should return:
#' #      Eta-squared adj Eta-squared         SE        LL        UL
#' # [1,]       0.241       0.2213707 0.06258283 0.1040229 0.3493431
#'  
#' 
#' @importFrom stats pf
#' @importFrom stats qnorm
#' @export
ci.etasqr <- function(alpha, etasqr, df1, df2) {
 alpha1 <- alpha/2
 alpha2 <- 1 - alpha1
 z0 <- qnorm(1 - alpha1)
 z <- qnorm(alpha2)
 F <- (etasqr/(1 - etasqr))*(df2/df1)
 adj <- 1 - (df2 + df1)*(1 - etasqr)/df2
 if (adj < 0) {adj = 0}
 ul0 <- 1 - exp(log(1 - etasqr) - z*sqrt(4*etasqr/(df2 - 1)))
 du <- ul0*(df1 + df2 + 1)/(1 - ul0)
 nc <- seq(0, du, by  = .001)
 p <- pf(F, df1, df2, nc)
 k1 <- which(min(abs(p - alpha2)) == abs(p - alpha2))[[1]]
 dl <- nc[k1]
 ll <- dl/(dl + df1 + df2 + 1)
 k2 <- which(min(abs(p - alpha1)) == abs(p - alpha1))[[1]]
 du <- nc[k2]
 ul <- du/(du + df1 + df2 + 1)
 if (ul == 0) {ul = ul0}
 se <- (ul - ll)/(2*z0)
 out <- t(c(etasqr, adj, se, ll, ul))
 colnames(out) <- c("Eta-squared", "adj Eta-squared", "SE", "LL", "UL")
 return(out)
}

# ci.2x2.mean.mixed ===========================================================
#' Computes tests and confidence intervals of effects in a 2x2 mixed design 
#' for means
#'
#'
#' @description
#' Computes confidence intervals and tests for the AB interaction effect, 
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
#' Computes confidence intervals and tests for the AB interaction effect, 
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
#' Computes tests and confidence intervals of effects in a 2x2 between-subjects 
#' design for means
#'
#'
#' @description
#' Computes confidence intervals and tests for the AB interaction effect, 
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


# ci.2x2.stdmean.bs ============================================================
#' Computes confidence intervals of standardized effects in a 2x2 
#' between-subjects design for means 
#'
#'
#' @description
#' Computes confidence intervals for standardized linear constrasts of means
#' (AB interaction, main effect of A, main efect of B, simple main effects
#' of A, and simple main effects of B) in a 2x2 between-subjects design with  
#' a quantitative response variable. Equality of population variances is not 
#' assumed. An unweigthed variance standardizer is used, which is the 
#' recommended standarizer when both factors are treatment factors.
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
#' * Estimate - bias adjusted estimate of standardized effect
#' * SE - standard error 
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @examples
#' y11 = c(14, 15, 11, 7, 16, 12, 15, 16, 10, 9)
#' y12 = c(18, 24, 14, 18, 22, 21, 16, 17, 14, 13)
#' y21 = c(16, 11, 10, 17, 13, 18, 12, 16, 6, 15)
#' y22 = c(18, 17, 11, 9, 9, 13, 18, 15, 14, 11)
#' ci.2x2.stdmean.bs(.05, y11, y12, y21, y22)
#'
#' # Should return:
#' #            Estimate        SE         LL         UL
#' # AB:      -1.4193502 0.6885238 -2.7992468 -0.1002829
#' # A:        0.4592015 0.3379520 -0.1933321  1.1314153
#' # B:       -0.7375055 0.3451209 -1.4297338 -0.0768846
#' # A at b1: -0.2504736 0.4640186 -1.1653006  0.6536189
#' # A at b2:  1.1688767 0.5001423  0.2136630  2.1741850
#' # B at a1: -1.4471806 0.4928386 -2.4441376 -0.5122457
#' # B at a2: -0.0278304 0.4820369 -0.9732017  0.9163482
#'
#'
#' @importFrom stats qnorm
#' @export
ci.2x2.stdmean.bs <- function(alpha, y11, y12, y21, y22) {
 z <- qnorm(1 - alpha/2)
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
 s <- sqrt((sd11^2 + sd12^2 + sd21^2 + sd22^2)/4)
 m <- c(m11, m12, m21, m22)
 sd <- c(sd11, sd12, sd21, sd22) 
 n <- c(n11, n12, n21, n22)
 var <- sd^2
 a <- 4
 df <- sum(n) - a
 adj <- 1 - 3/(4*df - 1)
 # AB 
 est1 <- (t(v1)%*%m)/s
 est1u <- adj*est1
 a1 <- est1^2/(2*a^2*s^4)
 a2 <- a1*sum((var^2/(n - 1)))
 a3 <- sum((v1^2*var/(n - 1)))/s^2
 se1 <- sqrt(a2 + a3)
 LL1 <- est1 - z*se1
 UL1 <- est1 + z*se1
 row1 <- c(est1u, se1, LL1, UL1)
# A 
 est2 <- (t(v2)%*%m)/s
 est2u <- adj*est2
 a1 <- est2^2/(2*a^2*s^4)
 a2 <- a1*sum((var^2/(n - 1)))
 a3 <- sum((v2^2*var/(n - 1)))/s^2
 se2 <- sqrt(a2 + a3)
 LL2 <- est2 - z*se2
 UL2 <- est2 + z*se2
 row2 <- c(est2u, se2, LL2, UL2)
# B 
 est3 <- (t(v3)%*%m)/s
 est3u <- adj*est3
 a1 <- est3^2/(2*a^2*s^4)
 a2 <- a1*sum((var^2/(n - 1)))
 a3 <- sum((v3^2*var/(n - 1)))/s^2
 se3 <- sqrt(a2 + a3)
 LL3 <- est3 - z*se3
 UL3 <- est3 + z*se3
 row3 <- c(est3u, se3, LL3, UL3)
# A at b1 
 est4 <- (t(v4)%*%m)/s
 est4u <- adj*est4
 a1 <- est4^2/(2*a^2*s^4)
 a2 <- a1*sum((var^2/(n - 1)))
 a3 <- sum((v4^2*var/(n - 1)))/s^2
 se4 <- sqrt(a2 + a3)
 LL4 <- est4 - z*se4
 UL4 <- est4 + z*se4
 row4 <- c(est4u, se4, LL4, UL4)
# A at b2 
 est5 <- (t(v5)%*%m)/s
 est5u <- adj*est5
 a1 <- est5^2/(2*a^2*s^4)
 a2 <- a1*sum((var^2/(n - 1)))
 a3 <- sum((v5^2*var/(n - 1)))/s^2
 se5 <- sqrt(a2 + a3)
 LL5 <- est5 - z*se5
 UL5 <- est5 + z*se5
 row5 <- c(est5u, se5, LL5, UL5)
# B at a1 
 est6 <- (t(v6)%*%m)/s
 est6u <- adj*est6
 a1 <- est6^2/(2*a^2*s^4)
 a2 <- a1*sum((var^2/(n - 1)))
 a3 <- sum((v6^2*var/(n - 1)))/s^2
 se6 <- sqrt(a2 + a3)
 LL6 <- est6 - z*se6
 UL6 <- est6 + z*se6
 row6 <- c(est6u, se6, LL6, UL6)
# B at a2 
 est7 <- (t(v7)%*%m)/s
 est7u <- adj*est7
 a1 <- est7^2/(2*a^2*s^4)
 a2 <- a1*sum((var^2/(n - 1)))
 a3 <- sum((v7^2*var/(n - 1)))/s^2
 se7 <- sqrt(a2 + a3)
 LL7 <- est7 - z*se7
 UL7 <- est7 + z*se7
 row7 <- c(est7u, se7, LL7, UL7)
 out <- rbind(row1, row2, row3, row4, row5, row6, row7)
 rownames(out) <- c("AB:", "A:", "B:", "A at b1:", "A at b2:", "B at a1:", "B at a2:")
 colnames(out) = c("Estimate", "SE", "LL", "UL")
 return(out)
}


# ci.2x2.median.bs =============================================================
#' Computes tests and confidence intervals of effects in a 2x2 betwen-subjects 
#' design for medians
#'
#'
#' @description
#' Computes confidence intervals for the AB interaction effect, main effect of
#' A, main efect of B, simple main effects of A, and simple main effects of B 
#' in a 2x2 between-subjects design with a quantitative response variable. The
#' effects are defined in terms of medians rather than means.
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
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @examples
#' y11 = c(14, 15, 11, 7, 16, 12, 15, 16, 10, 9)
#' y12 = c(18, 24, 14, 18, 22, 21, 16, 17, 14, 13)
#' y21 = c(16, 11, 10, 17, 13, 18, 12, 16, 6, 15)
#' y22 = c(18, 17, 11, 9, 9, 13, 18, 15, 14, 11)
#' ci.2x2.median.bs(.05, y11, y12, y21, y22)
#'
#' # Should return:
#' #          Estimate       SE         LL         UL
#' # AB:          -5.0 3.389735 -11.643758 1.64375833
#' # A:            1.5 1.694867  -1.821879 4.82187916
#' # B:           -2.0 1.694867  -5.321879 1.32187916
#' # A at b1:     -1.0 2.152661  -5.219138 3.21913797
#' # A at b2:      4.0 2.618464  -1.132095 9.13209504
#' # B at a1:     -4.5 2.311542  -9.030539 0.03053939
#' # B at a2:      0.5 2.479330  -4.359397 5.35939682
#'
#'
#' @importFrom stats qnorm
#' @importFrom stats pbinom
#' @importFrom stats median
#' @export
ci.2x2.median.bs <- function(alpha, y11, y12, y21, y22) {
 zcrit <- qnorm(1 - alpha/2)
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
 m11 <- median(y11)
 m12 <- median(y12)
 m21 <- median(y21)
 m22 <- median(y22)
 m <- c(m11, m12, m21, m22)
 n <- c(n11, n12, n21, n22)
 y11 <- sort(y11)
 y12 <- sort(y12)
 y21 <- sort(y21)
 y22 <- sort(y22)
 a <- round(n11/2 - sqrt(n11))
 if (a < 1) {a = 1}
 ll <- y11[a]
 ul <- y11[n11 - a + 1]
 p <- pbinom(a - 1, size = n11, prob = .5)
 z0 <- qnorm(1 - p)
 se11 <- (ul - ll)/(2*z0)
 a <- round(n12/2 - sqrt(n12))
 if (a < 1) {a = 1}
 ll <- y12[a]
 ul <- y12[n12 - a + 1]
 p <- pbinom(a - 1, size = n12, prob = .5)
 z0 <- qnorm(1 - p)
 se12 <- (ul - ll)/(2*z0)
 a <- round(n21/2 - sqrt(n21))
 if (a < 1) {a = 1}
 ll <- y21[a]
 ul <- y21[n21 - a + 1]
 p <- pbinom(a - 1, size = n21, prob = .5)
 z0 <- qnorm(1 - p)
 se21 <- (ul - ll)/(2*z0)
 a <- round(n22/2 - sqrt(n22))
 if (a < 1) {a = 1}
 ll <- y22[a]
 ul <- y22[n22 - a + 1]
 p <- pbinom(a - 1, size = n22, prob = .5)
 z0 <- qnorm(1 - p)
 se22 <- (ul - ll)/(2*z0)
 se <- c(se11, se12, se21, se22)
 # AB
 est1 <- t(v1)%*%m
 se1 <- sqrt(t(v1)%*%diag(se^2)%*%v1)
 LL1 <- est1 - zcrit*se1
 UL1 <- est1 + zcrit*se1
 row1 <- c(est1, se1, LL1, UL1)
 # A
 est2 <- t(v2)%*%m
 se2 <- sqrt(t(v2)%*%diag(se^2)%*%v2)
 LL2 <- est2 - zcrit*se2
 UL2 <- est2 + zcrit*se2
 row2 <- c(est2, se2, LL2, UL2)
 # B
 est3 <- t(v3)%*%m
 se3 <- sqrt(t(v3)%*%diag(se^2)%*%v3)
 LL3 <- est3 - zcrit*se3
 UL3 <- est3 + zcrit*se3
 row3 <- c(est3, se3, LL3, UL3)
 # A at b1
 est4 <- t(v4)%*%m
 se4 <- sqrt(t(v4)%*%diag(se^2)%*%v4)
 LL4 <- est4 - zcrit*se4
 UL4 <- est4 + zcrit*se4
 row4 <- c(est4, se4, LL4, UL4)
 # A at b2
 est5 <- t(v5)%*%m
 se5 <- sqrt(t(v5)%*%diag(se^2)%*%v5)
 LL5 <- est5 - zcrit*se5
 UL5 <- est5 + zcrit*se5
 row5 <- c(est5, se5, LL5, UL5)
 # B at a1
 est6 <- t(v6)%*%m
 se6 <- sqrt(t(v6)%*%diag(se^2)%*%v6)
 LL6 <- est6 - zcrit*se6
 UL6 <- est6 + zcrit*se6
 row6 <- c(est6, se6, LL6, UL6)
 # B at a2
 est7 <- t(v7)%*%m
 se7 <- sqrt(t(v7)%*%diag(se^2)%*%v7)
 LL7 <- est7 - zcrit*se7
 UL7 <- est7 + zcrit*se7
 row7 <- c(est7, se7, LL7, UL7)
 out <- rbind(row1, row2, row3, row4, row5, row6, row7)
 rownames(out) <- c("AB:", "A:", "B:", "A at b1:", "A at b2:", "B at a1:", "B at a2:")
 colnames(out) = c("Estimate", "SE", "LL", "UL")
 return(out)
}


#  =========================== Hypothesis tests ==============================
#  test.skew =================================================================
#' Computes p-value for test of skewness
#'
#'                          
#' @description
#' Computes a Monte Carlo p-value (250,000 replications) for the null 
#' hypothesis that the sample data come from a normal distribution. If the
#' p-value is small (e.g., less than .05) and the skewness estimate is 
#' positive, then the normality assumption can be rejected due to positive
#' skewness. If the p-value is small (e.g., less than .05) and the skewness
#' estimate is negative, then the normality assumption can be rejected due 
#' to negative skewness.
#'
#'
#' @param   y      vector of quantitative scores
#'  
#' 
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Skewness - estimate of skewness coefficient
#' * p - Monte Carlo two-sided p-value for test of zero skewness
#'
#'
#' @examples
#' y <- c(30, 20, 15, 10, 10, 60, 20, 25, 20, 30, 10, 5, 50, 40, 95)
#' test.skew(y)
#'
#' # Should return:
#' #      Skewness      p
#' # [1,]   1.5201 0.0067
#'
#'
#' @importFrom stats rnorm
#' @importFrom stats sd
#' @export
test.skew <- function(y) {
 rep <- 250000
 n = length(y)
 a <- sqrt((n - 1)/n)
 m <- mean(y)
 sd <- a*sd(y)
 z <- (y - m)/sd
 skew <- mean(z^3)
 y0 <- matrix(rnorm(rep*n), nrow = rep)
 m0 <- matrix(rowMeans(y0), nrow = rep, ncol = 1)
 sd0 <- a*sqrt(apply(y0, 1, var))
 sd0 <- matrix(sd0, nrow = rep, ncol = 1)
 j <- replicate(n, 1)
 z0 <- (y0 - m0%*%j)/(sd0%*%j)
 skew0 <- rowMeans(z0^3)
 c1 <- as.integer(skew0 > abs(skew))
 c2 <- as.integer(skew0 < -abs(skew))
 e1 <- sum(c1)/rep
 e2 <- sum(c2)/rep
 p <- e1 + e2
 out <- round(t(c(skew, p)), 4)
 colnames(out) <- c("Skewness", "p")
 return(out)
}


#  test.kurtosis =================================================================
#' Computes p-value for test of excess kurtosis
#'
#'                                 
#' @description
#' Computes a Monte Carlo p-value (250,000 replications) for the null 
#' hypothesis that the sample data come from a normal distribution. If the
#' p-value is small (e.g., less than .05) and excess kurtosis is positive,
#' then the normality assumption can be rejected due to leptokurtosis. If the
#' p-value is small (e.g., less than .05) and excess kurtosis is negative,
#' then the normality assumption can be rejected due to platykurtosis.
#'
#'
#' @param   y      vector of quantitative scores
#'  
#' 
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Kurtosis - estimate of kurtosis coefficient
#' * Excess - estimate of excess kurtosis (kurtosis - 3)
#' * p - Monte Carlo two-sided p-value for test of zero excess kurtosis
#'
#'
#' @examples
#' y <- c(30, 20, 15, 10, 10, 60, 20, 25, 20, 30, 10, 5, 50, 40, 95)
#' test.kurtosis(y)
#'
#' # Should return:
#' #      Kurtosis  Excess      p
#' # [1,]   4.8149  1.8149 0.0385 
#'
#'
#' @importFrom stats rnorm
#' @importFrom stats sd
#' @export
test.kurtosis <- function(y) {
 rep <- 250000
 n = length(y)
 a <- sqrt((n - 1)/n)
 m <- mean(y)
 sd <- a*sd(y)
 z <- (y - m)/sd
 kur <- mean(z^4)
 y0 <- matrix(rnorm(rep*n), nrow = rep)
 m0 <- matrix(rowMeans(y0), nrow = rep, ncol = 1)
 sd0 <- a*sqrt(apply(y0, 1, var))
 sd0 <- matrix(sd0, nrow = rep, ncol = 1)
 j <- replicate(n, 1)
 z0 <- (y0 - m0%*%j)/(sd0%*%j)
 kur0 <- rowMeans(z0^4)
 if (kur > 3) {c <- as.integer(kur0 > kur)} 
 if (kur < 3) {c <- as.integer(kur0 < kur)}
 p <- 2*sum(c)/rep
 if (p > .9999) {p = .9999}
 out <- round(t(c(kur, kur - 3, p)), 4)
 colnames(out) <- c("Kurtosis", "Excess", "p")
 return(out)
}


#  test.mono.mean.bs ==========================================================
#' Test of a monotonic trend in means for an ordered between-subjects factor
#' 
#'                     
#' @description
#' Computes simultaneous confidence intervals for all adjacent pairwise
#' comparisons of population means using estimated group means, estimated 
#' group standard deviations, and samples sizes as input. Equal variances are 
#' not assumed. A Satterthwaite adjustment to the degrees of freedom is used  
#' to improve the accuracy of the confidence intervals. If one or more lower
#' limits are greater than 0 and no upper limit is less than 0, then conclude
#' that the population means are monotoic decreasing. If one or more upper 
#' limits are less than 0 and no lower limits are greater than 0, then
#' conclude that the population means are monotoic increasing. Reject the 
#' hypothesis of a monotonic trend if any lower limit is greater than 0 and 
#' any upper limit is less than 0. 
#'
#'
#' @param  alpha   alpha level for simultaneous 1-alpha confidence
#' @param  m       vector of estimated group means
#' @param  sd      vector of estimated group standard deviations
#' @param  n       vector of sample sizes
#'
#'
#' @return 
#' Returns a matrix with the number of rows equal to the number
#' of adjacent pairwise comparisons. The columns are:
#' * Estimate - estimated mean difference
#' * SE - standard error
#' * LL - one-sided lower limit of the confidence interval
#' * UL - one-sided upper limit of the confidence interval
#'
#'
#' @examples
#' m <- c(12.86, 24.57, 36.29, 53.21)
#' sd <- c(13.185, 12.995, 14.773, 15.145)
#' n <- c(20, 20, 20, 20)
#' test.mono.mean.bs(.05, m, sd, n)
#'
#' # Should return:
#' #     Estimate       SE        LL         UL
#' # 1 2   -11.71 4.139530 -22.07803 -1.3419744
#' # 2 3   -11.72 4.399497 -22.74731 -0.6926939
#' # 3 4   -16.92 4.730817 -28.76921 -5.0707936
#'
#'
#' @importFrom stats qt
#' @export
test.mono.mean.bs <-function(alpha, m, sd, n) {
 a <- length(m)
 v <- sd^2
 m1 <- m[1: a - 1]
 m2 <- m[2: a]
 Estimate <- m1 - m2
 v1 <- v[1: a - 1]
 v2 <- v[2: a]
 n1 <- n[1: a - 1]
 n2 <- n[2: a]
 SE <- sqrt(v1/n1 + v2/n2)
 t <- Estimate/SE
 df <- SE^4/(v1^2/(n1^2*(n1 - 1)) + v2^2/(n2^2*(n2 - 1)))
 tcrit <- qt(1 - alpha/(2*(a - 1)), df)
 LL <- Estimate - tcrit*SE
 UL <- Estimate + tcrit*SE
 pair = cbind(seq(1, a - 1), seq(2, a))
 out <- cbind(pair, Estimate, SE, LL, UL)
 rownames(out) <- rep("", a - 1)
 return(out)
}


# ====================== Sample Size for Desired Precision ====================
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


#  size.ci.stdmean2 ===========================================================
#' Sample size for a 2-group standardized mean difference confidence interval
#'
#'
#' @description
#' Computes the sample size per group required to estimate two types of 
#' population standardized mean differences (unweighted standardizer and single
#' group standardizer) with desired confidence interval precision in a 2-group 
#' design. Set the standardized mean difference planning value to the largest 
#' value within a plausible range for a conservatively large sample size. Set
#' R = 1 for equal sample sizes.
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
#' Returns the required sample size per each group for each standardizer
#'
#'
#' @examples
#' size.ci.stdmean2(.05, .75, .5, 1)
#'
#' # Should return:
#' #                              n1  n2
#' # Unweighted standardizer:    132 132
#' # Single group standardizer:  141 141
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.ci.stdmean2 <- function(alpha, d, w, R) {
 z <- qnorm(1 - alpha/2)
 n11 <- ceiling((d^2*(1 + R)/(2*R) + 4*(1 + R)/R)*(z/w)^2)
 n12 <- ceiling(R*n11)
 n21 <- ceiling((2*d^2*(1 + R)/(2*R) + 4*(1 + R)/R)*(z/w)^2)
 n22 <- ceiling(R*n21)
 out1 <- t(c(n11, n12))
 out2 <- t(c(n21, n22))
 out <- rbind(out1, out2)
 colnames(out) <- c("n1", "n2")
 rownames(out) <- c("Unweighted standardizer:", "Single group standardizer:")
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
#' Sample size for a between-subjects standardized linear contrast of means
#' confidence interval
#'
#'                     
#' @description
#' Computes the sample size per group (assuming equal sample sizes) 
#' required to estimate two types of standardized linear contrasts of 
#' population means (unweighted average standardizer and single group
#' standardizer) with desired confidence interval precision in a 
#' between-subjects design. Set the standardized linear contrast of
#' means to the largest value within a plausible range for a conservatively
#' large sample size. 
#'
#'
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  d      planning value of standardized linear contrast of means 
#' @param  w      desired confidence interval width
#' @param  v      vector of between-subjects contrast coefficients 
#'
#'
#' @references
#' \insertRef{Bonett2009}{statpsych}
#'
#'
#' @return 
#' Returns the required sample size per group for each standardizer
#'
#'
#' @examples
#' v <- c(.5, .5, -.5, -.5)
#' size.ci.lc.stdmean.bs(.05, 1, .6, v)
#'
#' # Should return:
#' #                            Sample size per group
#' # Unweighted standardizer:                      49
#' # Single group standardizer:                    65
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.ci.lc.stdmean.bs <- function(alpha, d, w, v) {
 z <- qnorm(1 - alpha/2)
 a <- length(v)
 n1 <- ceiling((2*d^2/a + 4*(t(v)%*%v))*(z/w)^2)
 n2 <- ceiling((2*d^2 + 4*(t(v)%*%v))*(z/w)^2)
 out1 <- n1
 out2 <- n2
 out <- rbind(out1, out2)
 colnames(out) <- "Sample size per group"
 rownames(out) <- c("Unweighted standardizer:", "Single group standardizer:")
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
#' variances for the two measurements. Set the Pearson correlation planning  
#' value to the smallest value within a plausible range for a conservatively  
#' large sample size. Set the variance planning value to the largest value 
#' within a plausible range for a conservatively large sample size. 
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
#' Computes the sample size required to estimate two types of population 
#' standardized mean differences (unweighted standardizer and single group
#' standardizer) with desired confidence interval precision in a paired-samples
#' design. For a conservatively large sample size, set the standardized mean 
#' difference planning value to the largest value within a plausible range and 
#' set the Pearson correlation planning value to the smallest value within a 
#' plausible range.
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
#' Returns the required sample size for each standardizer 
#'
#'
#' @examples
#' size.ci.stdmean.ps(.05, 1, .65, .6)
#'
#' # Should return:
#' #                            Sample Size
#' # Unweighted standardizer:            46
#' # Single group standardizer:          52
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.ci.stdmean.ps <- function(alpha, d, cor, w) {
 z <- qnorm(1 - alpha/2)
 n1 <- ceiling((d^2*(1 + cor^2) + 8*(1 - cor))*(z/w)^2)
 n2 <- ceiling((2*d^2 + 8*(1 - cor))*(z/w)^2)
 out1 <- n1
 out2 <- n2
 out <- rbind(out1, out2)
 colnames(out) <- "Sample size"
 rownames(out) <- c("Unweighted standardizer:", "Single group standardizer:")
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
#' sample size. Set the Pearson correlation planning value to the 
#' smallest value within a plausible range for a conservatively 
#' large sample size. 
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


#  size.ci.lc.stdmean.ws ===================================================
#' Sample size for a within-subjects standardized linear contrast of means
#' confidence interval
#'
#'
#' @description
#' Computes the sample size required to estimate two types of standardized 
#' linear contrasts of population means (unweighted standardizer and single
#' level standardizer) with desired confidence interval precision in a 
#' within-subjects design. For a conservatively large sample size, set the 
#' standardized linear contrast of means planning value to the largest value
#' within a plausible range, and set the Pearson correlation planning value
#' to the smallest value within a plausible range.
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
#' Returns the required sample size for each standardizer
#'
#'
#' @examples
#' q <- c(.5, .5, -.5, -.5)
#' size.ci.lc.stdmean.ws(.05, 1, .7, .6, q)
#'
#' # Should return:
#' #                            Sample size
#' # Unweighted standardizer:            26
#' # Single level standardizer:          35
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.ci.lc.stdmean.ws <- function(alpha, d, cor, w, q) {
 z <- qnorm(1 - alpha/2)
 a <- length(q)
 n1 <- ceiling((2*d^2*(1 + (a - 1)*cor^2)/a + 4*(1 - cor)*(t(q)%*%q))*(z/w)^2)
 n2 <- ceiling((2*d^2 + 4*(1 - cor)*(t(q)%*%q))*(z/w)^2)
 out1 <- n1
 out2 <- n2
 out <- rbind(out1, out2)
 colnames(out) <- "Sample size"
 rownames(out) <- c("Unweighted standardizer:", "Single level standardizer:")
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
#' @param  r      number of measurements
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
size.ci.cronbach <- function(alpha, rel, r, w) {
 z <- qnorm(1 - alpha/2)
 n0 <- ceiling((8*r/(r - 1))*(1 - rel)^2*(z/w)^2 + 2)
 df1 = n0 - 1
 df2 = n0*(r - 1)
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


#  size.ci.etasqr =============================================================
#' Sample size for an eta-squared confidence interval 
#'
#'     
#' @description
#' Computes the sample size required to estimate an eta-squared coefficient
#' in a one-way ANOVA with desired confidence interval precision. Set the 
#' planning value of eta-squared to 1/3 for a conservatively large sample
#' size. 
#'
#'  
#' @param  alpha    alpha level for 1-alpha confidence
#' @param  etasqr   planning value of squared multiple correlation
#' @param  groups   number of groups
#' @param  w        desired confidence interval width
#'
#' 
#' @return 
#' Returns the total sample size requirement
#' 
#' 
#' @examples
#' size.ci.etasqr(.05, .333, 3, .2)
#'
#' # Should return:
#' #      Total sample size
#' # [1,]               189
#'  
#' 
#' @importFrom stats qnorm
#' @export                                                                  
size.ci.etasqr <- function(alpha, etasqr, groups, w) {
 alpha1 <- alpha/2
 alpha2 <- 1 - alpha1
 df1 <- groups - 1
 z <- qnorm(alpha2)
 n1 <- 16*etasqr*(1 - etasqr)*(z/w)^2 + groups + 1
 df2 <- n1 - groups
 ci <- ci.etasqr(alpha, etasqr, df1, df2)
 ll <- ci[1,4]                                  
 ul <- ci[1,5]                                  
 n2 <- ceiling(n1*((ul - ll)/w)^2)
 df2 <- n2 - groups
 ci <- ci.etasqr(alpha, etasqr, df1, df2)
 ll <- ci[1,4]                                  
 ul <- ci[1,5]
 n <- ceiling(n2*((ul - ll)/w)^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Total sample size"
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


# ======================== Sample Size for Desired Power ======================
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
#' Computes the sample size in each group required  to test a difference in 
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
#' with desired power in a paired-samples design. Set the Pearson correlation  
#' planning value to the smallest value within a plausible range, and set the
#' variance planning value to the largest value within a plausible range for a
#' conservatively large sample size.
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
#' alpha (e.g., use alpha = .10 rather alpha = .05). Set the Pearson correlation 
#' value to the smallest value within a plausible range, and set the variance
#' planning value to the largest value within a plausible range for a
#' conservatively large sample size.
#'
#'
#' @param   alpha  alpha level for hypothesis test 
#' @param   pow    desired power
#' @param   var    planning value of average variance of the two measurements
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
#' es > -h.  Set the Pearson correlation planning value to the smallest value
#' within a plausible range, and set the variance planning value to the largest
#' value within a plausible range for a conservatively large sample size.
#'
#'
#' @param   alpha  alpha level for hypothesis test 
#' @param   pow    desired power
#' @param   var    planning value of average variance of the two measurements
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
 es <- p - .5
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
 es <- p - .5
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
#' @param  r      number of measurements
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
size.test.cronbach <- function(alpha, pow, rel, r, h) {
 za <- qnorm(1 - alpha/2)
 zb <- qnorm(pow)
 e <- (1 - rel)/(1 - h)
 n <- ceiling((2*r/(r - 1))*(za + zb)^2/log(e)^2 + 2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 return(out)
}

# ======================= Power for Planned Sample Size =======================
#  power.mean1 ================================================================
#' Approximates the power of a one-sample t-test for a planned sample size
#'
#'
#' @description
#' Computes the approximate power of a one-sample t-test for a planned sample
#' size. For a conservatively low power approximation, set the variance 
#' planning value to the largest value within its plausible range, and set the 
#' effect size to a minimally interesting value.
#'
#'
#' @param  alpha  alpha level for hypothesis test 
#' @param  n      planned sample size
#' @param  var    planning value of response variable variance 
#' @param  es     planning value of mean minus hypothesized value
#'
#'
#' @return
#' Returns the approximate power of the test
#'
#'
#' @examples
#' power.mean1(.05, 15, 80.5, 7)
#'
#' # Should return:
#' #          Power
#' # [1,] 0.8021661
#'
#'
#' @importFrom stats qt
#' @importFrom stats pt
#' @export
power.mean1 <- function(alpha, n, var, es) { 
 df <- n - 1
 t <- qt(1 - alpha/2, df)
 z <- abs(es)/sqrt(var/n)
 pow <- 1 - pt(t, df, z)
 out <- matrix(pow, nrow = 1, ncol = 1)
 colnames(out) <- "Power"
 return(out)
}


#  power.mean2 ================================================================
#' Approximates the power of a two-sample t-test for planned sample sizes
#'
#'
#' @description
#' Computes the approximate power of a two-sample t-test for planned sample
#' sizes. For a conservatively low power approximation, set the variance 
#' planning values to the largest values within their plausible ranges, 
#' and set the effect size to a minimally interesting value. The within-group 
#' variances can be unequal across groups and a Satterthwaite degree of freedom 
#' adjustment is used to improve the accuracy of the power approximation.
#'
#'
#' @param  alpha  alpha level for hypothesis test 
#' @param  n1     planned sample size for group 1
#' @param  n2     planned sample size for group 2
#' @param  var1   planning value of within-group variance for group 1
#' @param  var2   planning value of within-group variance for group 2
#' @param  es     planning value of mean difference
#'
#'
#' @return
#' Returns the approximate power of the test
#'
#'
#' @examples
#' power.mean2(.05, 25, 25, 5.0, 6.0, 2)
#'
#' # Should return:
#' #          Power
#' # [1,] 0.8398413
#'
#'
#' @importFrom stats qt
#' @importFrom stats pt
#' @export
power.mean2 <- function(alpha, n1, n2, var1, var2, es) {
 df <- (var1/n1 + var2/n2)^2/(var1^2/(n1^2*(n1 - 1)) + var2^2/(n2^2*(n2 - 1)))
 t <- qt(1 - alpha/2, df)
 z <- abs(es)/sqrt(var1/n1 + var2/n2)
 pow <- 1 - pt(t, df, z)
 out <- matrix(pow, nrow = 1, ncol = 1)
 colnames(out) <- "Power"
 return(out)
}


#  power.lc.mean.bs ===========================================================
#' Approximates the power of a test for a linear contrast of means for planned 
#' sample sizes in a between-subjects design
#'
#'
#' @description
#' Computes the approximate power of a test for a linear contrast of population
#' mean for planned sample sizes in a between-subject design. The groups can be
#' the factor levels of a single factor design or the combinations of factors 
#' in a factorial design. For a conservatively low power approximation, set the
#' variance planning values to the largest values within their plausible ranges,
#' and set the effect size to a minimally interesting value. The within-group 
#' variances can be unequal across groups and a Satterthwaite degree of freedom 
#' adjustment is used to improve the accuracy of the power approximation.
#'
#'
#' @param  alpha  alpha level for hypothesis test 
#' @param  n      vector of planned sample sizes
#' @param  var    vector of within-group variance planning values
#' @param  es     planning value of linear contrast of means
#' @param  v      vector of contrast coefficients
#'
#'
#' @return
#' Returns the approximate power of the test
#'
#'
#' @examples
#' n <- c(20, 20, 20, 20)
#' var <- c(70, 70, 80, 80)
#' v <- c(.5, .5, -.5, -.5)
#' power.lc.mean.bs(.05, n, var, 5, v)
#'
#' # Should return:
#' #          Power
#' # [1,] 0.7221139
#'
#'
#' @importFrom stats qt
#' @importFrom stats pt
#' @export
power.lc.mean.bs <- function(alpha, n, var, es, v) {
 S <- diag(var)%*%(solve(diag(n)))
 se <- sqrt(t(v)%*%S%*%v)
 df <- (se^4)/sum(((v^4)*(var^2)/(n^2*(n - 1))))
 t <- qt(1 - alpha/2, df)
 z <- abs(es)/se
 pow <- 1 - pt(t, df, z)
 out <- matrix(pow, nrow = 1, ncol = 1)
 colnames(out) <- "Power"
 return(out)
}


#  power.mean.ps ==============================================================
#' Approximates the power of a paired-samples t-test for a planned sample size
#'
#'
#' @description
#' Computes the approximate power of a paired-samples t-test for a planned 
#' sample size. For a conservatively low power approximation, set the variance 
#' planning values to the largest values within their plausible ranges, set the
#' correlation planning value to the smallest value within its plausible range, 
#' and set the effect size to a minimally interesting value. The variances of
#' the two measurements can be unequal.
#'
#'
#' @param  alpha  alpha level for hypothesis test 
#' @param  n      planned sample size 
#' @param  var1   planning value of measurement 1 variance
#' @param  var2   planning value of measurement 2 variance
#' @param  es     planning value of mean difference
#' @param  cor    planning value of correlation between measurements
#'
#'
#' @return
#' Returns the approximate power of the test
#'
#'
#' @examples
#' power.mean.ps(.05, 20, 10.0, 12.0, 2, .7)
#'
#' # Should return:
#' #          Power
#' # [1,] 0.9074353
#'
#'
#' @importFrom stats qt
#' @importFrom stats pt
#' @export
power.mean.ps <- function(alpha, n, var1, var2, es, cor) {
 df <- n - 1
 t <- qt(1 - alpha/2, df)
 z <- abs(es)/sqrt((var1 + var2 - 2*cor*sqrt(var1*var2))/n)
 pow <- 1 - pt(t, df, z)
 out <- matrix(pow, nrow = 1, ncol = 1)
 colnames(out) <- "Power"
 return(out)
}


# ============================== Miscellaneous ===============================
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
#' @param  m      estimated mean
#' @param  sd     estimated standard deviation
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
#' Prediction interval for a difference of scores in a 2-group experiment
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
#' @param  m1     estaimted mean for group 1
#' @param  m2     estimated mean for group 1
#' @param  sd1    estimated standard deviation for group 1
#' @param  sd2    estimated standard deviation for group 2
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


#  pi.score.ps ================================================================ 
#' Prediction interval for difference of scores in a 2-level within-subjects 
#' experiment
#'
#'
#' @description
#' For a 2-level within-subjects experiment, this function computes a 
#' prediction interval for how the response variable score for one randomly 
#' selected person from the study population would differ under the two 
#' treatment conditions. 
#'
#'
#' @param  alpha  alpha level for 1-alpha confidence 
#' @param  m1     estimated mean from group 1
#' @param  m2     estimated mean from group 2
#' @param  sd1    estimated standard deviation from group 1
#' @param  sd2    estimated standard deviation from group 2
#' @param  cor    estimated correlation of paired scores
#' @param  n      sample size
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Predicted - predicted difference in scores
#' * df - degrees of freedom
#' * LL - lower limit of the prediction interval
#' * UL - upper limit of the prediction interval
#'
#'
#' @examples
#' pi.score.ps(.05, 265.1, 208.6, 23.51, 19.94, .814, 30)
#'
#' # Should return:
#' #      Predicted df       LL       UL
#' # [1,]      56.5 29 28.05936 84.94064
#'  
#' 
#' @importFrom stats qt
#' @export
pi.score.ps <- function(alpha, m1, m2, sd1, sd2, cor, n) {
 df <- n - 1
 tcrit <- qt(1 - alpha/2, df)
 est <- m1 - m2
 var <- sd1^2 + sd2^2 - 2*cor*sd1*sd2
 se <- sqrt(var + var/n)
 ll <- est - tcrit*se
 ul <- est + tcrit*se
 out <- t(c(est, df, ll, ul))
 colnames(out) <- c("Predicted", "df", "LL", "UL")
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
#' specified population mean and standard deviation. This functions is useful 
#' for generating hypothetical data for classroom demonstrations.
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
#' Computes an upper limit for a population variance using an estimated variance 
#' from a sample of size n in a prior study. The upper limit can be used as 
#' a variance planning value in sample size functions that require a variance 
#' planning value.
#'
#'
#' @param  alpha  alpha value for 1-alpha confidence (one-sided)
#' @param  var    estimated variance
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


#  pi.var.upper =============================================================== 
#' Upper prediction limit for an estimated variance
#'
#'                        
#' @description
#' Computes an approximate upper prediction limit for the estimated variance 
#' in a future study for a planned sample size. The prediction limit uses a 
#' variance estimate from a prior study. The upper variance prediction limit
#' is useful as a variance planning value for the sample size required to obtain
#' a confidence interval with desired width. This strategy for specifying a 
#' variance planning value is useful in applications where the population variance
#' in the prior study is assumed to be similar to the population variance in the
#' planned study. This variance planning value can be used to revise the planned 
#' sample size in the future study.
#'
#'
#' @param  alpha  alpha value for upper 1-alpha confidence 
#' @param  var    estimated variance from prior study
#' @param  n1     sample size used to estimate variance
#' @param  n2     planned sample size of future study
#'
#'
#' @return 
#' Returns an upper prediction estimate (UL) of an estimated variance in a future study
#'
#'
#' @examples
#' pi.var.upper(.2, 15, 40, 100)
#'
#' # Should return:
#' #           UL
#' # [1] 18.78522
#'  
#' 
#' @importFrom stats qnorm
#' @export
pi.var.upper <- function(alpha, var, n1, n2) {
 z <- qnorm(1 - alpha)
 ul <- exp(log(var) + z*sqrt(2/(n1 - 1) + 2/(n2 - 1)))
 out <- matrix(ul, nrow = 1, ncol = 1)
 colnames(out) <- c("UL")
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
 if (adj < 0) {adj = 0}
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
#' the estimated group means, estimated group standard deviations, and group
#' sample sizes.  
#'
#'
#' @param   m       vector of estimated group means
#' @param   sd      vector of estimated group standard deviations
#' @param   n       vector of group sample sizes
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


#  etasqr.gen.2way =============================================================== 
#' Generalized eta-squared estimates in a two-factor design
#'
#'
#' @description
#' Computes generalized eta-square estimates in a two-factor design where one
#' or both factors are classification factors. If both factors are treatment
#' factors, then partial eta-square estimates are typically recommended.
#' Use the estimates from this function in the etasqr.adj function to obtain
#' bias adjusted estimates.
#'
#'
#' @param  SSa    sum of squares for factor A
#' @param  SSb    sum of squares for factor B
#' @param  SSab   sum of squares for A x B interaction
#' @param  SSe    error (within) sum of squares
#'
#'
#' @return 
#' Returns a 3-row matrix. The columns are:
#' * A - estimate of eta-squared for factor A
#' * B - estimate of eta-squared for factor B
#' * AB - estimate of eta-squared for A x B interaction
#'
#'
#' @examples
#' etasqr.gen.2way(12.3, 15.6, 5.2, 7.9)
#'
#' # Should return:
#' #                                           A         B        AB
#' # A treatment, B classification:     0.300000 0.5435540 0.1811847
#' # A classification, B treatment:     0.484252 0.3804878 0.2047244
#' # A classification, B classiciation: 0.300000 0.3804878 0.1268293
#'  
#' 
#' @export
etasqr.gen.2way <- function(SSa, SSb, SSab, SSe) {
 etaA1 <- SSa/(SSa + SSb + SSab + SSe)
 etaB1 <- SSb/(SSb + SSab + SSe)
 etaAB1 <- SSab/(SSb + SSab + SSe)
 etaA2 <- SSa/(SSa + SSab + SSe)
 etaB2 <- SSb/(SSa + SSb + SSab + SSe)
 etaAB2 <- SSab/(SSa + SSab + SSe)
 etaA3 <- SSa/(SSa + SSb + SSab + SSe)
 etaB3 <- SSb/(SSa + SSb + SSab + SSe)
 etaAB3 <- SSab/(SSa + SSb + SSab + SSe)
 out1 <- t(c(etaA1, etaB1, etaAB1))
 out2 <- t(c(etaA2, etaB2, etaAB2))
 out3 <- t(c(etaA3, etaB3, etaAB3))
 out<- rbind(out1, out2, out3)
 colnames(out) <- c("A", "B", "AB")
 rownames1 <- c("A treatment, B classification:")
 rownames2 <- c("A classification, B treatment:")
 rownames3 <- c("A classification, B classiciation:")
 rownames(out) <- c(rownames1, rownames2, rownames3)
 return(out)
}


#  sim.ci.mean1 ===============================================================
#' Simulates confidence interval coverage probability for single mean
#'
#'                               
#' @description
#' Performs a computer simulation of the confidence interval performance for a 
#' single mean. Sample data can be generated from five different population 
#' distributions. All distributions are scaled to have standard deviations
#' of 1.0.
#'
#'
#' @param   alpha  alpha level for 1-alpha confidence
#' @param   n      sample size
#' @param   dist   type of distribution (1, 2, 3, 4,or 5)
#' * 1 = Gaussian (skewness = 0 and excess kurtosis = 0) 
#' * 2 = platykurtic (skewness = 0 and excess kurtosis = -1.2)
#' * 3 = leptokurtic (skewness = 0 and excess kurtsois = 6)
#' * 4 = moderate skew (skewness = 1 and excess kurtosis = 1.5)
#' * 5 = large skew (skewness = 2 and excess kurtosis = 6)
#' @param   rep    number of Monte Carlo samples
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
#' sim.ci.mean1(.05, 40, 4, 1000)
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
sim.ci.mean1 <- function(alpha, n, dist, rep) {
 tcrit <- qt(1 - alpha/2, n - 1)
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
#' Simulates confidence interval coverage probability for a 2-group mean 
#' difference
#'
#'                               
#' @description
#' Performs a computer simulation of separate variance and pooled variance 
#' confidence interval performance for a mean difference in a 2-group design.
#' Sample data within each group can be generated from five different 
#' population distributions. All distributions are scaled to have a standard 
#' deviation of 1.0 in group 1. 
#'
#'
#' @param   alpha     alpha level for 1-alpha confidence
#' @param   n1        sample size in group 1
#' @param   n2        sample size in group 2
#' @param   sd.ratio  ratio of population standard deviations (sd2/sd1)
#' @param   dist1     type of distribution for group 1 (1, 2, 3, 4, or 5)
#' @param   dist2     type of distribution for group 2 (1, 2, 3, 4, or 5)
#' * 1 = Gaussian (skewness = 0 and excess kurtosis = 0) 
#' * 2 = platykurtic (skewness = 0 and excess kurtosis = -1.2)
#' * 3 = leptokurtic (skewness = 0 and excess kurtsois = 6)
#' * 4 = moderate skew (skewness = 1 and excess kurtosis = 1.5)
#' * 5 = large skew (skewness = 2 and excess kurtosis = 6)
#' @param   rep       number of Monte Carlo samples
#'  
#' 
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Coverage - probability of confidence interval including population mean difference  
#' * Lower Error - probability of lower limit greater than population mean difference
#' * Upper Error - probability of upper limit less than population mean difference
#' * Ave CI Width - average confidence interval width
#'
#'
#' @examples
#' sim.ci.mean2(.05, 30, 25, 1.5, 4, 5, 1000)
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
sim.ci.mean2 <- function(alpha, n1, n2, sd.ratio, dist1, dist2, rep) {
 df1 <- n1 + n2 - 2
 tcrit1 <- qt(1 - alpha/2, df1)
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
#' Performs a computer simulation of confidence interval performance for a mean 
#' difference in a paired-samples design. Sample data within each level of the 
#' within-subjects factor can be generated from bivariate population
#' distributions with five different marginal distributions. All distributions
#' are scaled to have standard deviations of 1.0 at level 1. Bivariate random
#' data with specified marginal skewness and kurtosis are generated using the
#' unonr function in the mnonr package. 
#'
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
#' @param   rep       number of Monte Carlo samples
#'  
#' 
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Coverage - probability of confidence interval including population mean difference 
#' * Lower Error - probability of lower limit greater than population mean difference
#' * Upper Error - probability of upper limit less than population mean difference
#' * Ave CI Width - average confidence interval width
#'
#'
#' @examples
#' sim.ci.mean.ps(.05, 30, 1.5, .7, 4, 5, 1000)
#'
#' # Should return (within sampling error):
#' #      Coverage Lower Error Upper Error Ave CI Width
#' # [1,]  0.93815     0.05125      0.0106    0.7778518
#'
#'
#' @importFrom stats qt
#' @importFrom mnonr unonr
#' @export
sim.ci.mean.ps <- function(alpha, n, sd.ratio, cor, dist1, dist2, rep) {
 tcrit <- qt(1 - alpha/2, n - 1)
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
#' Performs a computer simulation of the confidence interval performance for 
#' a single median. Sample data can be generated from five different 
#' population distributions. All distributions are scaled to have standard 
#' deviations of 1.0.
#'
#'
#' @param   alpha  alpha level for 1-alpha confidence
#' @param   n      sample size
#' @param   dist   type of distribution (1, 2, 3, 4, or 5) 
#' * 1 = Gaussian (skewness = 0 and excess kurtosis = 0) 
#' * 2 = platykurtic (skewness = 0 and excess kurtosis = -1.2)
#' * 3 = leptokurtic (skewness = 0 and excess kurtsois = 6)
#' * 4 = moderate skew (skewness = 1 and excess kurtosis = 1.5)
#' * 5 = large skew (skewness = 2 and excess kurtosis = 6)
#' @param  rep     number of Monte Carlo samples
#'  
#' 
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Coverage - probability of confidence interval including population median  
#' * Lower Error - probability of lower limit greater than population mediann
#' * Upper Error - probability of upper limit less than population median
#' * Ave CI Width - average confidence interval width
#'
#'
#' @examples
#' sim.ci.median1(.05, 20, 5, 1000)
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
sim.ci.median1 <- function(alpha, n, dist, rep) {
 zcrit <- qnorm(1 - alpha/2)
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


#  sim.ci.median2 =============================================================
#' Simulates confidence interval coverage probability for a median difference
#' in a 2-group design
#'
#'                                          
#' @description
#' Performs a computer simulation of the confidence interval performance for a 
#' difference of medians in a 2-group design. Sample data for each group can
#' be generated from five different population distributions. All distributions
#' are scaled to have standard deviations of 1.0 in group 1.
#'
#' @param   alpha     alpha level for 1-alpha confidence
#' @param   n1        sample size for group 1
#' @param   n2        sample size for group 2
#' @param   sd.ratio  ratio of population standard deviations (sd2/sd1)
#' @param   dist1     type of distribution for group 1 (1, 2, 3, 4, or 5) 
#' @param   dist2     type of distribution for group 2 (1, 2, 3, 4, or 5) 
#' * 1 = Gaussian (skewness = 0 and excess kurtosis = 0) 
#' * 2 = platykurtic (skewness = 0 and excess kurtosis = -1.2)
#' * 3 = leptokurtic (skewness = 0 and excess kurtsois = 6)
#' * 4 = moderate skew (skewness = 1 and excess kurtosis = 1.5)
#' * 5 = large skew (skewness = 2 and excess kurtosis = 6)
#' @param   rep       number of Monte Carlo samples
#'  
#' 
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Coverage - Probability of confidence interval including population median difference  
#' * Lower Error - Probability of lower limit greater than population median difference
#' * Upper Error - Probability of upper limit less than population median difference
#' * Ave CI Width - Average confidence interval width
#'
#'
#' @examples
#' sim.ci.median2(.05, 20, 20, 2, 5, 4, 5000)
#'
#' # Should return (within sampling error):
#' #      Coverage Lower Error Upper Error Ave CI Width
#' # [1,]    0.952       0.027       0.021     2.368914
#'
#'
#' @importFrom stats qnorm
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom stats rt
#' @importFrom stats rgamma
#' @export
sim.ci.median2 <- function(alpha, n1, n2, sd.ratio, dist1, dist2, rep) {
 zcrit <- qnorm(1 - alpha/2)
 o1 <- round((n1/2 - sqrt(n1)))
 if (o1 < 1) {o1 = 1}
 o2 <- round((n2/2 - sqrt(n2)))
 if (o2 < 1) {o2 = 1}
 p1 <- pbinom(o1 - 1, size = n1, prob = .5)
 z01 <- qnorm(1 - p1)
 p2 <- pbinom(o2 - 1, size = n2, prob = .5)
 z02 <- qnorm(1 - p2)
 k <- 0; w <- 0; e1 <- 0; e2 <- 0
 repeat {
   k <- k + 1
   if (dist1 == 1) {
     y1 <- rnorm(n1)
     popmedian1 <- 0
   } else if (dist1 == 2) {
     y1 <- 3.464*runif(n1)
     popmedian1 <- 1.732
   } else if (dist1 == 3) {
     y1 <- .7745*rt(n1, 5)
     popmedian1 <- 0
   } else if (dist1 == 4) {
     y1 <- .5*rgamma(n1, 4)
     popmedian1 <- 1.837
   } else {
     y1 <- rgamma(n1, 1)
     popmedian1 <- 0.690
   }
   if (dist2 == 1) {
     y2 <- sd.ratio*rnorm(n2)
     popmedian2 <- 0
   } else if (dist2 == 2) {
     y2 <- sd.ratio*3.464*runif(n2)
     popmedian2 <- sd.ratio*1.732
   } else if (dist2 == 3) {
     y2 <- sd.ratio*.7745*rt(n2, 5)
     popmedian2 <- 0
   } else if (dist2 == 4) {
     y2 <- sd.ratio*.5*rgamma(n2, 4)
     popmedian2 <- sd.ratio*1.837
   } else {
     y2 <- sd.ratio*rgamma(n2, 1)
     popmedian2 <- sd.ratio*0.690
   }
   popdiff <- popmedian1 - popmedian2
   y1 <- sort(y1)
   LL1 <- y1[o1]
   UL1 <- y1[n1 - o1 + 1]
   median1 <- median(y1)
   se1 <- (UL1 - LL1)/(2*z01)
   y2 <- sort(y2)
   LL2 <- y2[o2]
   UL2 <- y2[n2 - o2 + 1]
   median2 <- median(y2)
   se2 <- (UL2 - LL2)/(2*z02)
   diff <- median1 - median2
   se <- sqrt(se1^2 + se2^2)
   ll <- diff - zcrit*se
   ul <- diff + zcrit*se
   w0 <- ul - ll
   c1 <- as.integer(ll > popdiff)
   c2 <- as.integer(ul < popdiff)
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


#  sim.ci.median.ps ===========================================================
#' Simulates confidence interval coverage probability for a median difference
#' in a paired-samples design
#'
#'                               
#' @description
#' Performs a computer simulation of confidence interval performance for a 
#' median difference in a paired-samples design. Sample data within each level
#' of the within-subjects factor can be generated from bivariate population
#' distributions with five different marginal distributions. All distributions
#' are scaled to have standard deviations of 1.0 at level 1. Bivariate random 
#' data with specified marginal skewness and kurtosis are generated using the
#' unonr function in the mnonr package. 
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
#' @param   rep       number of Monte Carlo samples
#'  
#' 
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Coverage - probability of confidence interval including population median difference  
#' * Lower Error - probability of lower limit greater than population median difference
#' * Upper Error - probability of upper limit less than population median difference
#' * Ave CI Width - average confidence interval width
#'
#'
#' @examples
#' sim.ci.median.ps(.05, 30, 1.5, .7, 4, 3, 2000)
#'
#' # Should return (within sampling error):
#' #      Coverage Lower Error Upper Error Ave CI Width
#' # [1,]    0.961       0.026       0.013    0.9435462
#'
#'
#' @importFrom stats qnorm
#' @importFrom mnonr unonr
#' @export
sim.ci.median.ps <- function(alpha, n, sd.ratio, cor, dist1, dist2, rep) {
 zcrit <- qnorm(1 - alpha/2)
 o <- round(n/2 - sqrt(n))
 if (o < 1) {o = 1}
 p <- pbinom(o - 1, size = n, prob = .5)
 z0 <- qnorm(1 - p)
 if (dist1 == 1) {
   skw1 <- 0; kur1 <- 0
   popmedian1 <- 0
 } else if (dist1 == 2) {
   skw1 <- 0; kur1 <- -1.2
   popmedian1 <- 0
 } else if (dist1 == 3) {
   skw1 <- 0; kur1 <- 6
   popmedian1 <- 0
 } else if (dist1 == 4) {
   skw1 <- .75; kur1 <- .86
   popmedian1 <- -.163
 } else {
   skw1 <- 1.41; kur1 <- 3
   popmedian1 <- -.313
 }
 if (dist2 == 1) {
   skw2 <- 0; kur2 <- 0
   popmedian2 <- 0
 } else if (dist2 == 2) {
   skw2 <- 0; kur2 <- -1.2
   popmedian2 <- 0
 } else if (dist2 == 3) {
   skw2 <- 0; kur2 <- 6
   popmedian2 <- 0
 } else if (dist2 == 4) {
   skw2 <- 1; kur2 <- 1.5
   popmedian2 <- -.163*sd.ratio
 } else {
   skw2 <- 2; kur2 <- 6
   popmedian2 <- -.313*sd.ratio
 }
 popdiff <- popmedian1 - popmedian2
 V <- matrix(c(1, cor*sd.ratio, cor*sd.ratio, sd.ratio^2), 2, 2)
 w <- 0; k <- 0; e1 <- 0; e2 <- 0
 repeat {
   k <- k + 1
   y <- unonr(n, c(0, 0), V, skewness = c(skw1, skw2), kurtosis = c(kur1, kur2))
   y1 <- y[, 1]
   y2 <- y[, 2]
   median1 <- median(y1)
   median2 <- median(y2)
   diff <- median1 - median2
   a1 <- (y1 < median1)
   a2 <- (y2 < median2)
   a3 <- a1 + a2
   a4 <- sum(a3 == 2)
   y1 <- sort(y1)
   y2 <- sort(y2)
   L1 <- y1[o]
   U1 <- y1[n - o + 1]
   se1 <- (U1 - L1)/(2*z0)
   L2 <- y2[o]
   U2 <- y2[n - o + 1]
   se2 <- (U2 - L2)/(2*z0)
   if (n/2 == trunc(n/2)) {
   p00 <- (sum(a4) + .25)/(n + 1)
   } else {
   p00 <- (sum(a4) + .25)/n 
   }
   cov <- (4*p00 - 1)*se1*se2
   se <- sqrt(se1^2 + se2^2 - 2*cov)
   ll <- diff - zcrit*se
   ul <- diff + zcrit*se
   w0 <- ul - ll
   c1 <- as.integer(ll > popdiff)
   c2 <- as.integer(ul < popdiff)
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


#  sim.ci.stdmean2 =============================================================
#' Simulates confidence interval coverage probability for a standardized mean
#' difference in a 2-group design
#'
#'                                      
#' @description
#' Performs a computer simulation of confidence interval performance for  
#' two types of standardized mean differences in a 2-group design. Sample 
#' data for each group can be generated from five different population 
#' distributions. All distributions are scaled to have standard deviations 
#' of 1.0 in group 1.
#'
#' @param   alpha     alpha level for 1-alpha confidence
#' @param   n1        sample size for group 1
#' @param   n2        sample size for group 2
#' @param   sd.ratio  ratio of population standard deviations (sd2/sd1)
#' @param   dist1     type of distribution for group 1 (1, 2, 3, 4, or 5) 
#' @param   dist2     type of distribution for group 2 (1, 2, 3, 4, or 5) 
#' * 1 = Gaussian (skewness = 0 and excess kurtosis = 0) 
#' * 2 = platykurtic (skewness = 0 and excess kurtosis = -1.2)
#' * 3 = leptokurtic (skewness = 0 and excess kurtsois = 6)
#' * 4 = moderate skew (skewness = 1 and excess kurtosis = 1.5)
#' * 5 = large skew (skewness = 2 and excess kurtosis = 6)
#' @param   d         population standardized mean difference 
#' @param   rep       number of Monte Carlo samples
#'  
#' 
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Coverage - Probability of confidence interval including population std mean difference  
#' * Lower Error - Probability of lower limit greater than population std mean difference
#' * Upper Error - Probability of upper limit less than population std mean difference
#' * Ave CI Width - Average confidence interval width
#'
#'
#' @examples
#' sim.ci.stdmean2(.05, 20, 20, 1.5, 3, 4, .75, 5000)
#'
#' # Should return (within sampling error):
#' #                         Coverage Lower Error Upper Error Ave CI Width   Ave Est
#' # Unweighted Standardizer   0.9058      0.0610      0.0332     1.342560 0.7838679
#' # Group 1 Standardizer      0.9450      0.0322      0.0228     1.827583 0.7862640
#'
#'
#' @importFrom stats qnorm
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom stats rt
#' @importFrom stats rgamma
#' @export
sim.ci.stdmean2 <- function(alpha, n1, n2, sd.ratio, dist1, dist2, d, rep) {
 zcrit <- qnorm(1 - alpha/2)
 df1 <- n1 - 1
 df2 <- n2 - 1
 df3 <- n1 + n2 - 2
 adj1 <- 1 - 3/(4*df3 - 1)
 adj2 <- 1 - 3/(4*df1 - 1)
 diff1 <- d*sqrt((1 + sd.ratio^2)/2)
 diff2 <- d
 k <- 0; w1 <- 0; w2 <- 0; e11 <- 0; e12 <- 0; e21 <- 0; e22 <- 0
 est1 <- 0; est2 <- 0;
 repeat {
   k <- k + 1
   if (dist1 == 1) {
     y1 <- rnorm(n1)  
   } else if (dist1 == 2) {
     y1 <- 3.464*runif(n1) - 1.732 
   } else if (dist1 == 3) {
     y1 <- .7745*rt(n1, 5) 
   } else if (dist1 == 4) {
     y1 <- .5*rgamma(n1, 4) - 2  
   } else {
     y1 <- rgamma(n1, 1) - 1 
   }
   if (dist2 == 1) {
     y0 <- sd.ratio*rnorm(n2)
   } else if (dist2 == 2) {
     y0 <- sd.ratio*3.464*runif(n2) - sd.ratio*1.734 
   } else if (dist2 == 3) {
     y0 <- sd.ratio*.7745*rt(n2, 5) 
   } else if (dist2 == 4) {
     y0 <- sd.ratio*.5*rgamma(n2, 4) - sd.ratio*2 
   } else {
     y0 <- sd.ratio*rgamma(n2, 1) - sd.ratio  
   }
   m1 <- mean(y1) + diff1
   m2 <- mean(y1) + diff2
   m0 <- mean(y0)
   v1 <- var(y1)
   v2 <- var(y0)
   s1 <- sqrt((v1 + v2)/2)
   d1 <- (m1 - m0)/s1
   se1 <- sqrt(d1^2*(v1^2/df1 + v2^2/df2)/(8*s1^4) + v1/(s1^2*df1) + v2/(s1^2*df2))
   ll1 <- d1 - zcrit*se1
   ul1 <- d1 + zcrit*se1
   est1 <- est1 + adj1*d1
   s2 <- sqrt(v1)
   d2 <- (m2 - m0)/s2
   se2 <- sqrt(d2^2/(2*df1) + 1/df1 + v2/(v1*df2))
   ll2 <- d2 - zcrit*se2
   ul2 <- d2 + zcrit*se2
   est2 <- est2 + adj2*d2
   w01 <- ul1 - ll1
   w02 <- ul2 - ll2
   w1 <- w1 + w01
   w2 <- w2 + w02
   c11 <- as.integer(ll1 > d)
   c21 <- as.integer(ul1 < d)
   c12 <- as.integer(ll2 > d)
   c22 <- as.integer(ul2 < d)
   e11 <- e11 + c11
   e21 <- e21 + c21
   e12 <- e12 + c12
   e22 <- e22 + c22
   if (k == rep) {break}
  }
  width1 <- w1/rep
  width2 <- w2/rep
  cov1 <- 1 - (e11 + e21)/rep
  cov2 <- 1 - (e12 + e22)/rep
  est1 = est1/rep
  est2 = est2/rep
  out1 <- t(c(cov1, e11/rep, e21/rep, width1, est1))
  out2 <- t(c(cov2, e12/rep, e22/rep, width2, est2))
  out <- rbind(out1, out2)
  colnames(out) <- c("Coverage", "Lower Error", "Upper Error", "Ave CI Width", "Ave Est")
  rownames(out) <- c("Unweighted Standardizer", "Group 1 Standardizer")
  return(out)
}


#  sim.ci.stdmean.ps ==========================================================
#' Simulates confidence interval coverage probability for a standardized mean
#' difference in a paired-samples design
#'
#'                                              
#' @description
#' Performs a computer simulation of confidence interval performance for  
#' two types of standardized mean differences in a paired-samples design. 
#' Sample data for each level of the within-subjects factor can be generated 
#' from five different population distributions. All distributions are scaled
#' to have standard deviations of 1.0 at level 1.
#'
#' @param   alpha     alpha level for 1-alpha confidence
#' @param   n         sample size 
#' @param   sd.ratio  ratio of population standard deviations (sd2/sd1)
#' @param   cor       correlation between paired measurements
#' @param   dist1     type of distribution at level 1 (1, 2, 3, 4, or 5) 
#' @param   dist2     type of distribution at level 2 (1, 2, 3, 4, or 5) 
#' * 1 = Gaussian (skewness = 0 and excess kurtosis = 0) 
#' * 2 = platykurtic (skewness = 0 and excess kurtosis = -1.2)
#' * 3 = leptokurtic (skewness = 0 and excess kurtsois = 6)
#' * 4 = moderate skew (skewness = 1 and excess kurtosis = 1.5)
#' * 5 = large skew (skewness = 2 and excess kurtosis = 6)
#' @param   d     population standardized mean difference 
#' @param   rep      number of Monte Carlo samples
#'  
#' 
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Coverage - Probability of confidence interval including population std mean difference
#' * Lower Error - Probability of lower limit greater than population std mean difference
#' * Upper Error - Probability of upper limit less than population std mean difference
#' * Ave CI Width - Average confidence interval width
#'
#'
#' @examples
#' sim.ci.stdmean.ps(.05, 20, 1.5, .8, 4, 4, .5, 2000)
#'
#' # Should return (within sampling error):
#' #                         Coverage Lower Error Upper Error Ave CI Width   Ave Est
#' # Unweighted Standardizer   0.9095      0.0555       0.035    0.7354865 0.5186796
#' # Level 1 Standardizer      0.9525      0.0255       0.022    0.9330036 0.5058198
#'
#'
#' @importFrom stats qnorm
#' @importFrom mnonr unonr
#' @export
sim.ci.stdmean.ps <- function(alpha, n, sd.ratio, cor, dist1, dist2, d, rep) {
 zcrit <- qnorm(1 - alpha/2)
 df <- n - 1
 adj1 <- sqrt((n - 2)/df)
 adj2 <- 1 - 3/(4*df - 1)
 diff1 <- d*sqrt((1 + sd.ratio^2)/2)
 diff2 <- d
 k <- 0; w1 <- 0; w2 <- 0; e11 <- 0; e12 <- 0; e21 <- 0; e22 <- 0
 est1 <- 0; est2 <- 0
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
 repeat {
   k <- k + 1
   y <- unonr(n, c(0, 0), V, skewness = c(skw1, skw2), kurtosis = c(kur1, kur2))
   y1 <- y[, 1] 
   y0 <- y[, 2]
   m1 <- mean(y1) + diff1
   m2 <- mean(y1) + diff2
   m0 <- mean(y0)
   v1 <- var(y1)
   v2 <- var(y0)
   s <- sqrt((v1 + v2)/2)
   cor <- cor(y1, y0)
   vd <- v1 + v2 - 2*cor*sqrt(v1*v2)
   d1 <- (m1 - m0)/s
   se1 <- sqrt(d1^2*(v1^2 + v2^2 + 2*cor^2*v1*v2)/(8*df*s^4) + vd/(df*s^2))
   ll1 <- d1 - zcrit*se1
   ul1 <- d1 + zcrit*se1
   est1 <- est1 + adj1*d1
   d2 <- (m2 - m0)/sqrt(v1)
   se2 <- sqrt(d2^2/(2*df) + vd/(df*v1))
   ll2 <- d2 - zcrit*se2
   ul2 <- d2 + zcrit*se2
   est2 <- est2 + adj2*d2
   w01 <- ul1 - ll1
   w02 <- ul2 - ll2
   c11 <- as.integer(ll1 > d)
   c21 <- as.integer(ul1 < d)
   c12 <- as.integer(ll2 > d)
   c22 <- as.integer(ul2 < d)
   e11 <- e11 + c11
   e21 <- e21 + c21
   e12 <- e12 + c12
   e22 <- e22 + c22
   w1 <- w1 + w01
   w2 <- w2 + w02
   if (k == rep) {break}
 }
 est1 = est1/rep
 est2 = est2/rep
 width1 <- w1/rep
 width2 <- w2/rep
 cov1 <- 1 - (e11 + e21)/rep
 cov2 <- 1 - (e12 + e22)/rep
 out1 <- t(c(cov1, e11/rep, e21/rep, width1, est1))
 out2 <- t(c(cov2, e12/rep, e22/rep, width2, est2))
 out <- rbind(out1, out2)
 colnames(out) <- c("Coverage", "Lower Error", "Upper Error", "Ave CI Width", "Ave Est")
 rownames(out) <- c("Unweighted Standardizer", "Level 1 Standardizer")
 return(out)
}


