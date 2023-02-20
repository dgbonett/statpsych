#  ======================= Confidence Intervals ==============================
#  ci.cor  ===================================================================
#' Confidence interval for a Pearson or partial correlation
#'
#'
#' @description
#' Computes a Fisher confidence interval for a population Pearson correlation  
#' or partial correlation with s control variables. Set s = 0 for a Pearson 
#' correlation. A bias adjustmentment is used to reduce the bias of the Fisher
#' transformed correlation. This function uses an estimated correlation as 
#' input. Use the cor.test function for raw data input.
#'
#'  
#' @param  alpha 	alpha level for 1-alpha confidence
#' @param  cor	  	estimated Pearson or partial correlation 
#' @param  s	  	number of control variables
#' @param  n	  	sample size
#'
#'
#' @references
#' \insertRef{Snedecor1980}{statpsych}
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimated correlation
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' ci.cor(.05, .536, 0, 50)
#'
#' # Should return:
#' #      Estimate        SE        LL        UL
#' # [1,]    0.536 0.1018149 0.2978573 0.7058914
#'  
#' 
#' @importFrom stats qnorm
#' @export
ci.cor <- function(alpha, cor, s, n) {
 z <- qnorm(1 - alpha/2)
 se <- sqrt((1 - cor^2)^2/(n - 1))
 se.z <- sqrt(1/((n - s - 3)))
 zr <- log((1 + cor)/(1 - cor))/2 - cor/(2*(n - 1))
 ll0 <- zr - z*se.z
 ul0 <- zr + z*se.z
 ll <- (exp(2*ll0) - 1)/(exp(2*ll0) + 1)
 ul <- (exp(2*ul0) - 1)/(exp(2*ul0) + 1)
 out <- t(c(cor, se, ll, ul))
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 return(out)
}


#  ci.spcor  =================================================================
#' Confidence interval for a semipartial correlation
#'
#'
#' @description
#' Computes a Fisher confidence interval for a population semipartial  
#' correlation. This function requires an (unadjusted) estimate of the 
#' squared multiple correlation in the full model that contains the
#' predictor variable of interest plus all control variables. This 
#' function computes a modified Aloe-Becker confidence interval that uses
#' n - 3 rather than n in the standard error and also uses a Fisher 
#' transformation of the semipartial correlation. 
#'
#'  
#' @param  alpha 	alpha level for 1-alpha confidence
#' @param  cor	  	estimated semipartial correlation 
#' @param  r2	  	estimated squared multiple correlation in full model
#' @param  n	  	sample size
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimated semipartial correlation
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @references
#' \insertRef{Aloe2012}{statpsych}
#'
#'
#' @examples
#' ci.spcor(.05, .582, .699, 20)
#'
#' # Should return:
#' #      Estimate        SE        LL        UL
#' # [1,]    0.582 0.1374298 0.2525662 0.7905182
#'  
#' 
#' @importFrom stats qnorm
#' @export
ci.spcor <- function(alpha, cor, r2, n) {
 z <- qnorm(1 - alpha/2)
 r0 <- r2 - cor^2
 zr <- log((1 + cor)/(1 - cor))/2
 a <- (r2^2 - 2*r2 + r0 - r0^2 + 1)/(1 - cor^2)^2
 se.z <- sqrt(a/(n - 3))
 se <- se.z*(1 - cor^2)
 ll0 <- zr - z*se.z
 ul0 <- zr + z*se.z
 ll <- (exp(2*ll0) - 1)/(exp(2*ll0) + 1)
 ul <- (exp(2*ul0) - 1)/(exp(2*ul0) + 1)
 out <- t(c(cor, se, ll, ul))
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 return(out)
}


#  ci.cor2 ====================================================================
#' Confidence interval for a 2-group Pearson correlation difference 
#'
#'
#' @description
#' Computes a confidence interval for a difference in population Pearson 
#' correlations in a 2-group design. A bias adjustmentment is used to 
#' reduce the bias of each Fisher transformed correlation.  
#'
#'  
#' @param  alpha	alpha level for 1-alpha confidence
#' @param  cor1	  	estimated Pearson correlation in group 1
#' @param  cor2	  	estimated Pearson correlation in group 2
#' @param  n1	  	sample size for group 1
#' @param  n2	  	sample size for group 2
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimated correlation difference
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @references
#' \insertRef{Zou2007}{statpsych}
#'
#'
#' @examples
#' ci.cor2(.05, .886, .802, 200, 200)
#'
#' # Should return:
#' #      Estimate         SE         LL        UL
#' # [1,]    0.084 0.02967934 0.02803246 0.1463609
#'  
#' 
#' @importFrom stats qnorm
#' @export
ci.cor2 <- function(alpha, cor1, cor2, n1, n2) {
 z <- qnorm(1 - alpha/2)
 diff <- cor1 - cor2
 se1 <- sqrt(1/(n1 - 3))
 zr1 <- log((1 + cor1)/(1 - cor1))/2 - cor1/(2*(n1 - 1))
 ll0 <- zr1 - z*se1
 ul0 <- zr1 + z*se1
 ll1 <- (exp(2*ll0) - 1)/(exp(2*ll0) + 1)
 ul1 <- (exp(2*ul0) - 1)/(exp(2*ul0) + 1)
 se2 <- sqrt(1/(n2 - 3))
 zr2 <- log((1 + cor2)/(1 - cor2))/2 - cor2/(2*(n2 - 1))
 ll0 <- zr2 - z*se2
 ul0 <- zr2 + z*se2
 ll2 <- (exp(2*ll0) - 1)/(exp(2*ll0) + 1)
 ul2 <- (exp(2*ul0) - 1)/(exp(2*ul0) + 1)
 ll <- diff - sqrt((cor1 - ll1)^2 + (ul2 - cor2)^2)
 ul <- diff + sqrt((ul1 - cor1)^2 + (cor2 - ll2)^2)
 se <- sqrt((1 - cor1^2)^2/(n1 - 3) + (1 - cor2^2)^2/((n2 - 3)))
 out <- t(c(diff, se, ll, ul))
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 return(out)
}


#  ci.cor.dep ================================================================
#' Confidence interval for a difference in dependent Pearson correlations
#'
#'
#' @description 
#' Computes a confidence interval for a difference in population Pearson 
#' correlations that are estimated from the same sample and have one 
#' variable in common. A bias adjustmentment is used to reduce the bias
#' of each Fisher transformed correlation. 
#'
#'  
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  cor1	  estimated Pearson correlation between y and x1
#' @param  cor2	  estimated Pearson correlation between y and x2
#' @param  cor12  estimated Pearson correlation between x1 and x2
#' @param  n	  sample size
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimated correlation difference
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @references
#' \insertRef{Zou2007}{statpsych}
#'
#'
#' @examples
#' ci.cor.dep(.05, .396, .179, .088, 166)
#'
#' # Should return:
#' #      Estimate         LL      UL
#' # [1,]    0.217 0.01323072 0.415802
#'  
#' 
#' @importFrom stats qnorm
#' @export 
ci.cor.dep <- function(alpha, cor1, cor2, cor12, n) {
 z <- qnorm(1 - alpha/2)
 diff <- cor1 - cor2
 se1 <- sqrt(1/((n - 3)))
 zr1 <- log((1 + cor1)/(1 - cor1))/2 - cor1/(2*(n - 1))
 ll0 <- zr1 - z*se1
 ul0 <- zr1 + z*se1
 ll1 <- (exp(2*ll0) - 1)/(exp(2*ll0) + 1)
 ul1 <- (exp(2*ul0) - 1)/(exp(2*ul0) + 1)
 se2 <- sqrt(1/((n - 3)))
 zr2 <- log((1 + cor2)/(1 - cor2))/2 - cor2/(2*(n - 1))
 ll0 <- zr2 - z*se2
 ul0 <- zr2 + z*se2
 ll2 <- (exp(2*ll0) - 1)/(exp(2*ll0) + 1)
 ul2 <- (exp(2*ul0) - 1)/(exp(2*ul0) + 1)
 b <- (1 - cor1^2)*(1 - cor2^2)
 c <- ((cor12 - cor1*cor2/2)*(1 - cor1^2 - cor2^2 - cor12^2) + cor12^3)/b
 d1 <- 2*c*(cor1 - ll1)*(cor2 - ul2)
 d2 <- 2*c*(cor1 - ul1)*(cor2 - ll2)
 ll <- diff - sqrt((cor1 - ll1)^2 + (ul2 - cor2)^2 - d1)
 ul <- diff + sqrt((ul1 - cor1)^2 + (cor2 - ll2)^2 - d2)
 out <- t(c(diff, ll, ul))
 colnames(out) <- c("Estimate", "LL", "UL")
 return(out)
}


#  ci.cor2.gen ================================================================
#' Confidence interval for a 2-group correlation difference
#'
#'
#' @description 
#' Computes a 100(1 - alpha)% confidence interval for a difference in 
#' population correlations in a 2-group design. The correlations can be 
#' Pearson, Spearman, partial, semipartial, or point-biserial correlations. 
#' The function requires 100(1 - alpha)% confidence intervals for each 
#' correlation as input. An approximate standard error is recovered from 
#' the confidence interval. 
#'
#'  
#' @param  cor1  estimated correlation for group 1 
#' @param  ll1   lower limit for group 1 correlation
#' @param  ul1   upper limit for group 1 correlation
#' @param  cor2  estimated correlation for group 2
#' @param  ll2   lower limit for group 2 correlation
#' @param  ul2   upper limit for group 2 correlation
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimated correlation difference
#' * SE - recovered standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @references
#' \insertRef{Zou2007}{statpsych}
#'
#'
#' @examples
#' ci.cor2.gen(.4, .35, .47, .2, .1, .32)
#'
#' # Should return:
#' #      Estimate   LL        UL
#' # [1,]      0.2 0.07 0.3220656
#'  
#' 
#' @importFrom stats qnorm
#' @export 
ci.cor2.gen <- function(cor1, ll1, ul1, cor2, ll2, ul2) {
 diff <- cor1 - cor2
 ll <- diff - sqrt((cor1 - ll1)^2 + (ul2 - cor2)^2)
 ul <- diff + sqrt((ul1 - cor1)^2 + (cor2 - ll2)^2)
 se <- (ul - ll)/(2*z)
 out <- t(c(diff, se, ll, ul))
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 return(out)
}


#  ci.pbcor ==================================================================
#' Confidence interval for a point-biserial correlation
#'
#'
#' @description 
#' Computes confidence intervals for two types of population point-biserial 
#' correlations. One type uses a weighted average of the group variances 
#' and is appropriate for nonexperimental designs with simple random sampling
#' (rather than stratified random sampling). The other type uses an unweighted 
#' average of the group variances and is appropriate for experimental designs.
#' Equality of variances is not assumed for either type. 
#'
#'  
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  m1     estimated mean for group 1
#' @param  m2     estimated mean for group 2
#' @param  sd1    estimated standard deviation for group 1
#' @param  sd2    estimated standard deviation for group 2
#' @param  n1     sample size for group 1
#' @param  n2	  sample size for group 2
#'
#'
#' @references
#' \insertRef{Bonett2020a}{statpsych}
#'
#'
#' @return 
#' Returns a 2-row matrix. The columns are:
#' * Estimate - estimated point-biserial correlation
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' ci.pbcor(.05, 28.32, 21.48, 3.81, 3.09, 40, 40)
#'
#' # Should return:
#' #              Estimate        LL        UL
#' # Weighted:   0.7065799 0.5885458 0.7854471
#' # Unweighted: 0.7020871 0.5808366 0.7828948
#'  
#' 
#' @importFrom stats qnorm
#' @export 
ci.pbcor <- function(alpha, m1, m2, sd1, sd2, n1, n2) {
 z <- qnorm(1 - alpha/2)
 s <- sqrt((sd1^2 + sd2^2)/2)
 d1 <- (m1 - m2)/sqrt(((n1 - 1)*sd1^2 + (n2 - 1)*sd2^2)/(n1 + n2 - 2))
 d2 <- (m1 - m2)/s
 n <- n1 + n2
 p <- n1/n
 k <- (n1 + n2 - 2)/(n*p*(1 - p))
 cor1 <- d1/sqrt(d1^2 + k)
 cor2 <- d2/sqrt(d2^2 + 4)
 sed1 <- sqrt(d1^2*(1/n1 + 1/n2)/8 + (sd1^2/n1 + sd2^2/n2)/s^2)
 sed2 <- sqrt(d2^2*(sd1^4/(n1-1) + sd2^4/(n2-1))/(8*s^4) + (sd1^2/(n1-1) + sd2^2/(n2-1))/s^2) 
 lld1 <- d1 - z*sed1
 uld1 <- d1 + z*sed1
 lld2 <- d2 - z*sed2
 uld2 <- d2 + z*sed2
 llc1 <- lld1/sqrt(lld1^2 + k)
 ulc1 <- uld1/sqrt(uld1^2 + k)
 llc2 <- lld2/sqrt(lld2^2 + 4)
 ulc2 <- uld2/sqrt(uld2^2 + 4)
 out1 <- t(c(cor1, llc1, ulc1))
 out2 <- t(c(cor2, llc2, ulc2))
 out <- rbind(out1, out2)
 colnames(out) <- c("Estimate", "LL", "UL")
 rownames(out) <- c("Weighted:", "Unweighted:")
 return(out)
}


#  ci.spear ==================================================================
#' Confidence interval for a Spearman correlation
#'
#'
#' @description
#' Computes a Fisher confidence interval for a population Spearman correlation.  
#'
#'  
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  y	  vector of y scores 
#' @param  x	  vector of x scores (paired with y)
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimated correlation
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @references
#' \insertRef{Bonett2000}{statpsych}
#'
#'
#' @examples
#' y <- c(21, 4, 9, 12, 35, 18, 10, 22, 24, 1, 6, 8, 13, 16, 19)
#' x <- c(67, 28, 30, 28, 52, 40, 25, 37, 44, 10, 14, 20, 28, 40, 51)
#' ci.spear(.05, y, x)
#'
#' # Should return:
#' #       Estimate         SE        LL        UL
#' # [1,] 0.8699639 0.08241326 0.5840951 0.9638297
#'  
#' 
#' @importFrom stats qnorm
#' @export
ci.spear <- function(alpha, y, x) {
 z <- qnorm(1 - alpha/2)
 n <- length(y)
 yr <- rank(y)
 xr <- rank(x)
 cor <- cor(yr, xr)
 se <- sqrt((1 + cor^2/2)*(1 - cor^2)^2/(n - 3))
 z.cor <- log((1 + cor)/(1 - cor))/2
 ll0 <- z.cor - z*se/(1 - cor^2)
 ul0 <- z.cor + z*se/(1 - cor^2)
 ll <- (exp(2*ll0) - 1)/(exp(2*ll0) + 1)
 ul <- (exp(2*ul0) - 1)/(exp(2*ul0) + 1)
 out <- cbind(cor, se, ll, ul)
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 return (out)
}
	
	
#  ci.spear2 ====================================================================
#' Confidence interval for a 2-group Spearman correlation difference 
#'
#'
#' @description
#' Computes a confidence interval for a difference in population Spearman 
#' correlations in a 2-group design. 
#'
#'  
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  cor1	  estimated Spearman correlation for group 1
#' @param  cor2	  estimated Spearman correlation for group 2
#' @param  n1	  sample size for group 1
#' @param  n2	  sample size for group 2
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimated correlation difference
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @references
#' \insertRef{Bonett2000}{statpsych}
#'
#'
#' \insertRef{Zou2007}{statpsych}
#'
#'
#' @examples
#' ci.spear2(.05, .54, .48, 180, 200)
#'
#' # Should return:
#' #      Estimate         SE         LL        UL
#' # [1,]     0.06 0.08124926 -0.1003977 0.2185085     
#'  
#' 
#' @importFrom stats qnorm
#' @export
ci.spear2 <- function(alpha, cor1, cor2, n1, n2) {
 z <- qnorm(1 - alpha/2)
 se1 <- sqrt((1 + cor1^2/2)*(1 - cor1^2)^2/(n1 - 3))
 se2 <- sqrt((1 + cor2^2/2)*(1 - cor2^2)^2/(n2 - 3))
 z.cor1 <- log((1 + cor1)/(1 - cor1))/2
 z.cor2 <- log((1 + cor2)/(1 - cor2))/2
 ll01 <- z.cor1 - z*se1/(1 - cor1^2)
 ul01 <- z.cor1 + z*se1/(1 - cor1^2)
 ll1 <- (exp(2*ll01) - 1)/(exp(2*ll01) + 1)
 ul1 <- (exp(2*ul01) - 1)/(exp(2*ul01) + 1)
 ll02 <- z.cor2 - z*se2/(1 - cor2^2)
 ul02 <- z.cor2 + z*se2/(1 - cor2^2)
 ll2 <- (exp(2*ll02) - 1)/(exp(2*ll02) + 1)
 ul2 <- (exp(2*ul02) - 1)/(exp(2*ul02) + 1)
 cor.dif <- cor1 - cor2
 ll.dif <- cor.dif - sqrt((cor1 - ll1)^2 + (ul2 - cor2)^2)
 ul.dif <- cor.dif + sqrt((ul1 - cor1)^2 + (cor2 - ll2)^2)
 se <- sqrt(se1^2 + se2^2)
 out <- cbind(cor.dif, se, ll.dif, ul.dif)
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 return (out)
}
	

#  ci.mape  ===================================================================
#' Confidence interval for a mean absolute prediction error
#'
#'
#' @description
#' Computes a confidence interval for a population mean absolute prediction
#' error (MAPE) in a general linear model. The MAPE is a more robust 
#' alternative to the residual standard deviation. This function requires a
#' vector of estimated residuals from a general linear model. This confidence
#' interval does not assume zero excess kurtosis but does assume symmetry of
#' the population prediction errors.
#'
#'  
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  res    vector of residuals 
#' @param  s	  number of predictor variables in model
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimated mean absolute prediction error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' res <- c(-2.70, -2.69, -1.32, 1.02, 1.23, -1.46, 2.21, -2.10, 2.56,
#'       -3.02, -1.55, 1.46, 4.02, 2.34)
#' ci.mape(.05, res, 1)
#'
#' # Should return:
#' #       Estimate       LL       UL
#' # [1,]    2.3744 1.751678 3.218499
#'  
#' 
#' @importFrom stats qt
#' @importFrom stats sd
#' @export
ci.mape <- function(alpha, res, s) {
 n <- length(res)
 df <- n - s - 1
 t <- qt(1 - alpha/2, df)
 c <- n/(n - (s + 2)/2)
 mape <- mean(abs(res))
 kur <- (sd(res)/mape)^2
 se <- sqrt((kur - 1)/df)
 ll <- exp(log(c*mape) - t*se)
 ul <- exp(log(c*mape) + t*se)
 out <- t(c(c*mape, ll, ul))
 colnames(out) <- c("Estimate", "LL", "UL")
 return(out)
}


#  ci.condslope  ==============================================================
#' Confidence interval for conditional (simple) slopes in a linear model
#'
#'
#' @description
#' Computes confidence intervals and test statistics for population 
#' conditional slopes (simple slopes) in a general linear model that
#' includes a predictor variable that is the product of a moderator 
#' variable and a predictor variable. Conditional slopes are computed  
#' at specified low and high values of the moderator variable. 
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
#' @param  dfe    error degrees of freedom 
#'
#'
#' @return 
#' Returns a 2-row matrix. The columns are:
#' * Estimate - estimated conditional slope
#' * t - t test statistic
#' * p - p-value
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' ci.condslope(.05, .132, .154, .031, .021, .015, 5.2, 10.6, 122)
#'
#' # Should return:
#' #                   Estimate        SE        t  df           p 
#' # At low moderator    0.9328 0.4109570 2.269824 122 0.024973618 
#' # At high moderator   1.7644 0.6070517 2.906507 122 0.004342076 
#' #                           LL       UL
#' # At low moderator   0.1192696 1.746330
#' # At high moderator  0.5626805 2.966119
#'  
#' 
#' @importFrom stats qt
#' @export                  
ci.condslope <- function(alpha, b1, b2, se1, se2, cov, lo, hi, dfe) {
 t <- qt(1 - alpha/2, dfe)
 slope.lo <- b1 + b2*lo
 slope.hi <- b1 + b2*hi
 se.lo <- sqrt(se1^2 + se2^2*lo^2 + 2*lo*cov)
 se.hi <- sqrt(se1^2 + se2^2*hi^2 + 2*hi*cov)
 t.lo <- slope.lo/se.lo
 t.hi <- slope.hi/se.hi
 p.lo <- 2*(1 - pt(abs(t.lo), dfe))
 p.hi <- 2*(1 - pt(abs(t.hi), dfe))
 ll.lo <- slope.lo - t*se.lo
 ul.lo <- slope.lo + t*se.lo
 ll.hi <- slope.hi - t*se.hi
 ul.hi <- slope.hi + t*se.hi
 out1 <- t(c(slope.lo, se.lo, t.lo, dfe, p.lo, ll.lo, ul.lo))
 out2 <- t(c(slope.hi, se.hi, t.hi, dfe, p.hi, ll.hi, ul.hi))
 out <- rbind(out1, out2)
 colnames(out) <- c("Estimate", "SE", "t", "df", "p", "LL", "UL")
 rownames(out) <- c("At low moderator", "At high moderator")
 return(out)
}


#  ci.lc.reg  ==============================================================
#' Confidence interval for a linear contrast of regression coefficients in
#' multiple group regression model
#'
#'  
#' @description
#' Compute a confidence interval and test statistic for a linear contrast
#' of a population regression coefficients (y-intercept or slope) across
#' groups in a multiple group regression model. Equality of error variances
#' across groups is not assumed. A Satterthwaite adjustment to the degrees 
#' of freedom is used to improve the accuracy of the confidence interval. 
#'
#'  
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  est    vector of parameter estimates
#' @param  se     vector of standard errors
#' @param  n      vector of group sample sizes
#' @param  s      number of predictor variables for each within-group model
#' @param  v      vector of contrast coefficients
#'
#' 
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimated linear contrast
#' * SE - standard error
#' * t - t test statistic
#' * df - degrees of freedom
#' * p - p-value
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' est <- c(1.74, 1.83, 0.482)
#' se <- c(.483, .421, .395)
#' n <- c(40, 40, 40)
#' v <- c(.5, .5, -1)
#' ci.lc.reg(.05, est, se, n, 4, v)
#'
#' # Should return:
#' #      Estimate        SE        t      df          p        LL       UL
#' # [1,]    1.303 0.5085838 2.562016 78.8197 0.01231256 0.2906532 2.315347
#'  
#' 
#' @importFrom stats qt
#' @export                  
ci.lc.reg <- function(alpha, est, se, n, s, v) {
 con <- t(v)%*%est 
 se.con <- sqrt(sum(v^2*se^2))
 t <- con/se.con
 df <- (sum(v^2*se^2)^2)/sum((v^4*se^4)/(n - s - 1))
 p <- 2*(1 - pt(abs(t), df))
 tcrit <- qt(1 - alpha/2, df)
 ll <- con - tcrit*se.con
 ul <- con + tcrit*se.con
 out <- t(c(con, se.con, t, df, p, ll, ul))
 colnames(out) <- c("Estimate", "SE", "t", "df", "p", "LL", "UL")
 return(out)
}


#  ci.fisher ================================================================= 
#' Fisher confidence interval 
#'
#' 
#' @description
#' Computes a Fisher confidence interval for any type of correlation or any
#' measure of association that has a -1 to 1 range.
#'
#'
#' @param   alpha  alpha value for 1-alpha confidence
#' @param   cor    estimated correlation or association coefficient 
#' @param   se     standard error of estimate
#'
#'
#' @return
#' Returns a 1-row matrix containing the lower and upper confidence limits.
#'
#'
#' @examples
#' ci.fisher(.05, .641, .052)
#'
#' # Should return:
#' #             LL        UL
#' # [1,] 0.5276396 0.7319293
#'
#'
#' @importFrom stats qnorm
#' @export
ci.fisher <- function(alpha, cor, se) {
 z <- qnorm(1 - alpha/2)
 zr <- log((1 + cor)/(1 - cor))/2
 ll0 <- zr - z*se/(1 - cor^2)
 ul0 <- zr + z*se/(1 - cor^2)
 ll <- (exp(2*ll0) - 1)/(exp(2*ll0) + 1)
 ul <- (exp(2*ul0) - 1)/(exp(2*ul0) + 1)
 out <- t(c(ll, ul))
 colnames(out) <- c("LL", "UL")
 return(out)
}


#  ci.indirect ===============================================================
#' Confidence interval for an indirect effect
#'
#'
#' @description
#' Computes a Monte Carlo confidence interval (500,000 trials) for a population
#' unstandardized indirect effect in a path model and a Sobel standard error. 
#' This function is not recommended for a standardized indirect effect unless 
#' the standardized slope estimates for both paths are less than about .3 in 
#' absolute value. The Monte Carlo method is general in that the slope estimates
#' and standard errors do not need to be OLS estimates with homoscedastic 
#' standard errors. For example, LAD slope estimates and their standard errors,
#' OLS slope estimates and heteroscedastic standard errors, distribution-free 
#' Theil-Sen slope estimates with McKean-Schrader standard errors, or 
#' standardized slopes with robust standard errors also could be used.
#'
#'  
#' @param  alpha  alpha level for 1-alpha confidence  
#' @param  b1     unstandardized slope estimate for first path
#' @param  b2     unstandardized slope estimate for second path
#' @param  se1    standard error for b1
#' @param  se2    standard error for b2
#'
#' 
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimated indirect effect
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' ci.indirect (.05, 2.48, 1.92, .586, .379)
#'
#' # Should return (within sampling error):
#' #      Estimate       SE       LL       UL
#' # [1,]   4.7616 1.625282 2.178812 7.972262
#'  
#' 
#' @importFrom stats rnorm
#' @export  
ci.indirect <- function(alpha, b1, b2, se1, se2) {
 k <- 500000
 c <- round(k*alpha/2)
 y1 <- rnorm(k, b1, se1)
 y2 <- rnorm(k, b2, se2)
 y <- sort(y1*y2)
 se <- sqrt((b1*se2)^2 + (b2*se1)^2)
 ll <- y[c]
 ul <- y[k - c + 1]
 out <- t(c(b1*b2, se, ll, ul))
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 return(out)
}


#  ci.lc.gml ==================================================================
#' Confidence interval for a linear contrast of general linear model parameters
#'
#'                                  
#' @description
#' Computes the estimate, standard error, and confidence interval for a linear
#' contrast of parameters in a general linear model using coef(object) and
#' vcov(object) where "object" is a fitted model object from the lm function.
#'
#'
#' @param   alpha  alpha for 1 - alpha confidence
#' @param   n      sample size
#' @param   b      vector of parameter estimates from coef(object)
#' @param   V      covariance matrix of parameter estimates from vcov(object)
#' @param   q      vector of coefficients
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimate of linear function 
#' * SE - standard error
#' * t - t test statistic 
#' * df - degrees of freedom
#' * p - p-value 
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval 
#'
#'
#' @examples
#' y <- c(43, 62, 49, 60, 36, 79, 55, 42, 67, 50)
#' x1 <- c(3, 6, 4, 6, 2, 7, 4, 2, 7, 5)
#' x2 <- c(4, 6, 3, 7, 1, 9, 3, 3, 8, 4)
#' out <- lm(y ~ x1 + x2)
#' b <- coef(out)
#' V <- vcov(out)
#' n <- length(y)
#' q <- c(0, .5, .5)
#' b
#' ci.lc.glm(.05, n, b, V, q)
#'
#' #  Should return:
#' # (Intercept)          x1          x2 
#' #   26.891111    3.648889    2.213333 
#' # > ci.lc.glm(.05, n, b, V, q)
#' #      Estimate        SE       t df           p       LL       UL
#' # [1,] 2.931111 0.4462518 6.56829  7 0.000313428 1.875893 3.986329
#'
#'
#' @importFrom stats qt
#' @importFrom stats pt
#' @export
ci.lc.glm <-function(alpha, n, b, V, q) {
 df <- n - length(b)
 tcrit <- qt(1 - alpha/2, df)
 est <- t(q)%*%b
 se <- sqrt(t(q)%*%V%*%q)
 t <- est/se
 p <- 2*(1 - pt(abs(t), df))
 ll <- est - tcrit*se
 ul <- est + tcrit*se
 out <- t(c(est, se, t, df, p, ll, ul))
 colnames(out) <- c("Estimate", "SE", "t", "df", "p", "LL", "UL")
 return(out)
}


#  ci.lc.gen.bs ===============================================================
#' Confidence interval for a linear contrast of parameters in a between-subjects
#' design
#'
#'                                              
#' @description
#' Computes the estimate, standard error, and approximate confidence interval 
#' for a linear contrast of any type of parameter (e.g., quartile, ordinal 
#' regression slope, path coefficient, G-index) where each parameter value has
#' been estimated from a different sample. The parameter vaues are assumed to 
#' be of the same type (e.g., all unstandardized path coefficients) and their 
#' sampling distributions are assumed to be approximately normal.
#'
#'
#' @param  alpha   alpha level for simultaneous 1-alpha confidence
#' @param  est     vector of parameter estimates
#' @param  se      vector of standard errors
#' @param  v       vector of contrast coefficients
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimate of linear contrast
#' * SE - standard error of linear contrast
#' * LL - lower limit of confidence interval
#' * UL - upper limit of confidence interval
#'
#'
#' @examples
#' est <- c(3.86, 4.57, 2.29, 2.88)
#' se <- c(0.185, 0.365, 0.275, 0.148)
#' v <- c(.5, .5, -.5, -.5)
#' ci.lc.gen.bs(.05, est, se, v)
#'
#' # Should return:
#' #      Estimate        SE       LL       UL
#' # [1,]     1.63 0.2573806 1.125543 2.134457
#'
#'
#' @importFrom stats qnorm
#' @export
ci.lc.gen.bs <- function(alpha, est, se, v) {
 est.lc <- t(v)%*%est
 se.lc <- sqrt(t(v)%*%diag(se^2)%*%v)
 zcrit <- qnorm(1 - alpha/2)
 ll <- est.lc - zcrit*se.lc
 ul <- est.lc + zcrit*se.lc
 out <- t(c(est.lc, se.lc, ll, ul))
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 return(out)
}


#  ci.rsqr ===================================================================
#' Confidence interval for squared multiple correlation
#'                             
#' @description
#' Computes an approximate confidence interval for a population squared 
#' multiple correlation in a linear model with random predictor variables.  
#' This function uses the scaled central F approximation method. An
#' approximate standard error is recovered from the confidence interval.
#'
#'
#' @param  alpha    alpha value for 1-alpha confidence
#' @param  r2       estimated unadjusted squared multiple correlation
#' @param  s        number of predictor variables
#' @param  n        sample size
#'
#' @references
#' \insertRef{Helland1987}{statpsych}
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * R-squared - estimate of unadjusted R-squared 
#' * adj R-squared - bias adjusted R-squared estimate
#' * SE - recovered standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @examples
#' ci.rsqr(.05, .241, 3, 116)
#'
#' # Should return:
#' #      R-squared adj R-squared         SE         LL        UL  
#' # [1,]     0.241     0.2206696 0.06752263 0.09819599 0.3628798
#'  
#' 
#' @importFrom stats qf
#' @export
ci.rsqr <- function(alpha, r2, s, n) {
 alpha1 <- alpha/2
 alpha2 <- 1 - alpha1
 z0 <- qnorm(1 - alpha1)
 dfe <- n - s - 1
 adj <- 1 - (n - 1)*(1 - r2)/dfe
 if (adj < 0) {adj = 0}
 b1 <- r2/(1 - r2)
 b2 <- adj/(1 - adj)
 v1 <- ((n - 1)*b1 + s)^2/((n - 1)*b1*(b1 + 2) + s)
 v2 <- ((n - 1)*b2 + s)^2/((n - 1)*b2*(b2 + 2) + s)
 F1 <- qf(alpha1, v1, dfe)
 F2 <- qf(alpha2, v2, dfe)
 ll <- (dfe*r2 - (1 - r2)*s*F2)/(dfe*(r2 + (1 - r2)*F2))
 ul <- (dfe*r2 - (1 - r2)*s*F1)/(dfe*(r2 + (1 - r2)*F1))
 if (ll < 0) {ll = 0}
 if (ul < 0) {ul = 0}
 i <- 1
 while (i < 30) {
   i <- i + 1
   b1 <- ul/(1 - ul)
   b2 <- ll/(1 - ll)
   v1 <- ((n - 1)*b1 + s)^2/((n - 1)*b1*(b1 + 2) + s)
   v2 <- ((n - 1)*b2 + s)^2/((n - 1)*b2*(b2 + 2) + s)
   F1 <- qf(alpha1, v1, dfe)
   F2 <- qf(alpha2, v2, dfe)
   ll <- (dfe*r2 - (1 - r2)*s*F2)/(dfe*(r2 + (1 - r2)*F2))
   ul <- (dfe*r2 - (1 - r2)*s*F1)/(dfe*(r2 + (1 - r2)*F1))
   if (ll < 0) {ll = 0}
   if (ul < 0) {ul = 0}
 }
 se <- (ul - ll)/(2*z0)
 out <- t(c(r2, adj, se, ll, ul))
 colnames(out) <- c("R-squared", "adj R-squared", "SE", "LL", "UL")
 return(out)
}


#  ci.theil ===================================================================
#' Theil-Sen estimate and confidence interval for slope
#'
#'                                 
#' @description
#' Computes a Theil-Sen estimate and distribution-free confidence interval 
#' for the slope of a simple linear regression model. An approximate 
#' standard error is recovered from the confidence interval.
#'
#'
#' @param   alpha   alpha level for 1-alpha confidence
#' @param   y       vector of response variable scores
#' @param   x       vector of predictor variable scores (paired with y)
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - Theil-Sen estimate of population slope
#' * SE - recovered standard error
#' * LL - lower limit of confidence interval
#' * UL - upper limit of confidence interval
#'
#'
#' @references
#' \insertRef{Hollander1999}{statpsych}
#'
#'
#' @examples
#' y <- c(21, 4, 9, 12, 35, 18, 10, 22, 24, 1, 6, 8, 13, 16, 19)
#' x <- c(67, 28, 30, 28, 52, 40, 25, 37, 44, 10, 14, 20, 28, 40, 51)
#' ci.theil(.05, y, x)
#'
#' # Should return:
#' #      Estimate        SE        LL   UL
#' # [1,]      0.5 0.1085927 0.3243243 0.75
#'
#'
#' @importFrom stats qnorm
#' @export
ci.theil <- function(alpha, y, x) {
  z <- qnorm(1 - alpha/2)
  n = length(x)
  x.p <- t(combn(x,2))
  y.p <- t(combn(y,2))
  y.d <- y.p[,1] - y.p[,2]
  x.d <- x.p[,1] - x.p[,2]
  s = which(x.d != 0, arr.ind = T)
  x.diff <- x.d[s]
  y.diff <- y.d[s]
  k <- length(x.diff)
  c = z*sqrt(k*(2*n + 5)/9) 
  o1 <- floor((k - c)/2)
  if (o1 < 1) {o1 = 1}
  o2 <- ceiling((k + c)/2 + 1)
  if (o2 > k) {o2 = k}
  b <- y.diff/x.diff
  b <- sort(b)
  est <- median(b)
  ll <- b[o1]
  ul <- b[o2]
  se <- (ul - ll)/(2*z)
  out <- t(c(est, se, ll, ul))
  colnames(out) = c("Estimate", "SE", "LL", "UL")
  return(out)
}


#  =================== Sample Size for Desire Precision =======================
#  size.ci.slope ==============================================================
#' Sample size for a slope confidence interval
#'
#'
#' @description
#' Computes the total sample size required to estimate a slope with desired 
#' confidence interval precision in a between-subjects design with a 
#' quantitative factor. In an experimental design, the total sample size 
#' would be allocated to the levels of the quantitative factor and it might
#' be necessary to increase the total sample size to achieve equal sample
#' sizes. Set the error variance planning value to the largest value within 
#' a plausible range for a conservatively large sample size.
#'
#'  
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  evar   planning value of within group (error) variance
#' @param  x      vector of x values of the quantitative factor
#' @param  w      desired confidence interval width
#'
#' 
#' @return 
#' Returns the required total sample size
#' 
#' 
#' @examples
#' x <- c(2, 5, 8)
#' size.ci.slope(.05, 31.1, x, 1)
#'
#' # Should return:
#' #      Sample size
#' # [1,]          83
#'  
#' 
#' @importFrom stats qnorm
#' @export  
size.ci.slope <- function(alpha, evar, x, w) {
 m <- mean(x)
 xvar <- sum((x - m)^2)/length(x)
 z <- qnorm(1 - alpha/2)
 n <- ceiling(4*(evar/xvar)*(z/w)^2 + 1 + z^2/2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 return(out)
}


#  size.ci.cor ==============================================================
#' Sample size for a Pearson or partial correlation confidence interval 
#'
#'
#' @description
#' Computes the sample size required to estimate a Pearson or
#' partial correlation with desired confidence interval precision. 
#' Set s = 0 for a Pearson correlation. Set the correlation planning value
#' to the smallest value within a plausible range for a conservatively 
#' large sample size.
#'
#'  
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  cor    planning value of correlation
#' @param  s      number of control variables 
#' @param  w      desired confidence interval width
#'
#' 
#' @references
#' \insertRef{Bonett2000}{statpsych}
#'
#'
#' @return 
#' Returns the required sample size
#' 
#' 
#' @examples
#' size.ci.cor(.05, .362, 0, .25)
#'
#' # Should return:
#' #      Sample size
#' # [1,]         188
#'  
#' 
#' @importFrom stats qnorm
#' @export  
size.ci.cor <- function(alpha, cor, s, w) {
 z <- qnorm(1 - alpha/2)
 n1 <- ceiling(4*(1 - cor^2)^2*(z/w)^2 + s + 3)
 zr <- log((1 + cor)/(1 - cor))/2
 se <- sqrt(1/(n1 - s - 3))
 ll0 <- zr - z*se
 ul0 <- zr + z*se
 ll <- (exp(2*ll0) - 1)/(exp(2*ll0) + 1)
 ul <- (exp(2*ul0) - 1)/(exp(2*ul0) + 1)
 n <- ceiling((n1 - s - 3)*((ul - ll)/w)^2 + 3 + s)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 return(out)
}


#  size.ci.rsqr ==============================================================
#' Sample size for a squared multiple correlation confidence interval
#'
#'
#' @description
#' Computes the sample size required to estimate a squared multiple correlation
#' in a random-x regression model with desired confidence interval precision.
#' Set the planning value of the squared multiple correlation to 1/3 for a 
#' conservatively large sample size. This function uses an approximation to
#' the standard error of the squared multiple correlation.
#'
#'  
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  r2     planning value of squared multiple correlation
#' @param  s      number of predictor variables in model
#' @param  w      desired confidence interval width
#'
#' 
#' @return 
#' Returns the required sample size
#' 
#' 
#' @examples
#' size.ci.rsqr(.05, .333, 2, .2)
#'
#' # Should return:
#' #      Sample size
#' # [1,]         232
#'  
#' 
#' @importFrom stats qnorm
#' @export  
size.ci.rsqr <- function(alpha, r2, s, w) {
 z <- qnorm(1 - alpha/2)
 n <- ceiling(16*(r2*(1 - r2)^2)*(z/w)^2 + s + 2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 return(out)
}


#  size.ci.condmean ==========================================================
#' Sample size for a conditional mean confidence interval  
#'
#'
#' @description
#' Computes the total sample size required to estimate a conditional mean of
#' y at x = x* in a fixed-x simple linear regression model with desired 
#' confidence interval precision. In an experimental design, the total sample
#' size would be allocated to the levels of the quantitative factor and it
#' might be necessary to increase the total sample size to achieve equal 
#' sample sizes. Set the error variance planning value to the largest value 
#' within a plausible range for a conservatively large sample size.
#'
#'  
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  evar   planning value of within group (error) variance
#' @param  xvar   variance of fixed predictor variables 
#' @param  diff   difference between x* and mean of x
#' @param  w      desired confidence interval width
#'
#' 
#' @return 
#' Returns the required total sample size
#' 
#' 
#' @examples
#' size.ci.condmean(.05, 120, 125, 15, 5)
#'
#' # Should return:
#' #      Total sample size
#' # [1,]               210
#'  
#' 
#' @importFrom stats qnorm
#' @export  
size.ci.condmean <- function(alpha, evar, xvar, diff, w) {
 z <- qnorm(1 - alpha/2)
 n <- ceiling(4*(evar*(1 + diff^2/xvar))*(z/w)^2 + 1 + z^2/2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Total sample size"
 return(out)
}


#  size.ci.lc.ancova =========================================================
#' Sample size for a linear contrast confidence interval in an ANCOVA  
#'
#'
#' @description
#' Computes the sample size for each group (assuming equal sample sizes) 
#' required to estimate a linear contrast of means in an ANCOVA model with 
#' desired confidence interval precision. In a nonexperimental design, the 
#' sample size is affected by the magnitude of covariate mean differences 
#' across groups. The covariate mean differences can be approximated by 
#' specifying the largest standardized covariate mean difference across all 
#' pairwise group differences and for all covariates. In an experiment, this
#' standardized mean difference should be set to 0. Set the error variance 
#' planning value to the largest value within a plausible range for a 
#' conservatively large sample size.
#'
#'  
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  evar   planning value of within group (error) variance
#' @param  s      number of covariates 
#' @param  d      largest standardized mean difference for all covariates
#' @param  w      desired confidence interval width
#' @param  v      vector of between-subjects contrast coefficients
#'
#' 
#' @return 
#' Returns the required sample size per group
#' 
#' 
#' @examples
#' v <- c(1, -1)
#' size.ci.lc.ancova(.05, 1.37, 1, 0, 1.5, v)
#'
#' # Should return:
#' #      Sample size per group
#' # [1,]                    21
#'  
#' 
#' @importFrom stats qnorm
#' @export  
size.ci.lc.ancova <- function(alpha, evar, s, d, w, v) {
 z <- qnorm(1 - alpha/2)
 k <- length(v)
 n <- ceiling(4*evar*(1 + d^2/4)*(t(v)%*%v)*(z/w)^2 + s + z^2/(2*k))
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size per group"
 return(out)
}


# ======================= Sample Size for Desired Power =======================
#  size.test.slope ============================================================
#' Sample size for a test of a slope
#'
#'
#' @description
#' Computes the total sample size required to test a population slope with  
#' desired power in a between-subjects design with a quantitative factor. 
#' In an experimental design, the total sample size would be allocated to the
#' levels of the quantitative factor and it might be necessary to use a larger
#' total sample size to achieve equal sample sizes. Set the error variance 
#' planning value to the largest value within a plausible range for a 
#' conservatively large sample size.
#'
#'  
#' @param  alpha   alpha level for hypothesis test
#' @param  pow     desired power
#' @param  evar    planning value of within-group (error) variance
#' @param  x       vector of x values of the quantitative factor
#' @param  slope   planning value of slope
#' @param  h       hypothesized value of slope  
#'
#' 
#' @return 
#' Returns the required total sample size
#' 
#' 
#' @examples
#' x <- c(2, 5, 8)
#' size.test.slope(.05, .9, 31.1, x, .75, 0)
#'
#' # Should return:
#' #      Sample size
#' # [1,]         100
#'  
#' 
#' @importFrom stats qnorm
#' @export  
size.test.slope <- function(alpha, pow, evar, x, slope, h) {
 za <- qnorm(1 - alpha/2)
 zb <- qnorm(pow)
 m <- mean(x)
 xvar <- sum((x - m)^2)/length(x)
 es <- slope - h
 n <- ceiling((evar/xvar)*(za + zb)^2/es^2 + 1 + za^2/2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 return(out)
}


#  size.test.cor ==============================================================
#' Sample size for a test of a Pearson or partial correlation 
#'
#'
#' @description
#' Computes the sample size required to test a Pearson or a partial correlation 
#' with desired power. Set s = 0 for a Pearson correlation. 
#'
#'  
#' @param  alpha   alpha level for hypothesis test
#' @param  pow     desired power
#' @param  cor     planning value of correlation
#' @param  s       number of control variables
#' @param  h       hypothesized value of correlation  
#'
#' 
#' @return 
#' Returns the required sample size
#' 
#' 
#' @examples
#' size.test.cor(.05, .9, .45, 0, 0)
#'
#' # Should return:
#' #      Sample size
#' # [1,]          48
#'  
#' 
#' @importFrom stats qnorm
#' @export  
size.test.cor <- function(alpha, pow, cor, s, h) {
 za <- qnorm(1 - alpha/2)
 zb <- qnorm(pow)
 zr <- log((1 + cor)/(1 - cor))/2
 zo <- log((1 + h)/(1 - h))/2
 es <- zr - zo
 n <- ceiling((za + zb)^2/es^2 + s + 3)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 return(out)
}


#  size.interval.cor =========================================================
#' Sample size for a finite interval test of a Pearson or partial correlation 
#'
#'
#' @description
#' Computes the sample size required to perform a finite interval test for a 
#' Pearson or a partial correlation with desired power. Set s = 0 for a 
#' Pearson correlation. The correlation planning value must be a value within 
#' the hypothesized finite interval. 
#'
#'  
#' @param  alpha   alpha level for hypothesis test
#' @param  pow     desired power
#' @param  cor     planning value of correlation
#' @param  s       number of control variables
#' @param  h       upper limit of hypothesized interval  
#'
#' 
#' @return 
#' Returns the required sample size
#' 
#' 
#' @examples
#' size.interval.cor(.05, .8, .1, 0, .25)
#'
#' # Should return:
#' #      Sample size
#' # [1,]         360
#'  
#' 
#' @importFrom stats qnorm
#' @export  
size.interval.cor <- function(alpha, pow, cor, s, h) {
 za <- qnorm(1 - alpha)
 zb <- qnorm(1 - (1 - pow)/2)
 zr <- log((1 + cor)/(1 - cor))/2
 zh <- log((1 + h)/(1 - h))/2
 es <- zh - zr
 n <- ceiling((za + zb)^2/es^2 + s + 3)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 return(out)
}


#  size.test.lc.ancova ==========================================================
#' Sample size for a mean linear contrast test in an ANCOVA 
#'
#'
#' @description
#' Computes the sample size for each group (assuming equal sample sizes) required
#' to test a linear contrast of means in an ANCOVA model with desired power. In a 
#' nonexperimental design, the sample size is affected by the magnitude of 
#' covariate mean differences across groups. The covariate mean differences can be 
#' approximated by specifying the largest standardized covariate mean difference 
#' across all pairwise comparisons and for all covariates. In an experiment, this 
#' standardized mean difference is set to 0. Set the error variance planning 
#' value to the largest value within a plausible range for a conservatively 
#' large sample size.
#'
#'  
#' @param  alpha   alpha level for hypothesis test
#' @param  pow     desired power
#' @param  evar    planning value of within-group (error) variance
#' @param  es      planning value of linear contrast
#' @param  s       number of covariates 
#' @param  d       largest standardized mean difference for all covariates
#' @param  v       vector of between-subjects contrast coefficients
#'
#' 
#' @return 
#' Returns the required sample size per group
#' 
#' 
#' @examples
#' v <- c(.5, .5, -1)
#' size.test.lc.ancova(.05, .9, 1.37, .7, 1, 0, v)
#'
#' # Should return:
#' #      Sample size per group
#' # [1,]                    47
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.test.lc.ancova <- function(alpha, pow, evar, es, s, d, v) {
 za <- qnorm(1 - alpha/2)
 zb <- qnorm(pow)
 k <- length(v)
 n <- ceiling((evar*(1 + d^2/4)*t(v)%*%v)*(za + zb)^2/es^2 + s + za^2/k)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size per group"
 return(out)
}


# ============================= Miscellaneous =================================
#  slope.contrast =============================================================
#' Contrast coefficients for the slope of a quantitative factor
#'
#'
#' @description
#' Computes the contrast coefficients to estimate the slope of a line in a
#' single factor design with a quantitative factor.
#'
#'
#' @param   x      vector of numeric factor levels
#'
#'
#' @return 
#' Returns the vector of contrast coefficients
#'
#'
#' @examples
#' x <- c(25, 50, 75, 100)
#' slope.contrast(x)
#'
#' # Should return:
#' #      Coefficient
#' # [1,]      -0.012
#' # [2,]      -0.004
#' # [3,]       0.004
#' # [4,]       0.012
#'  
#' 
#' @export
slope.contrast <- function(x) {
 a <- length(x)
 m <- matrix(1, a, 1)*mean(x)
 ss <- sum((x - m)^2)
 coef <- (x - m)/ss
 out <- matrix(coef, nrow = a, ncol = 1)
 colnames(out) = "Coefficient"
 return(out)
}


#  random.yx =================================================================
#' Generates random bivariate scores 
#'
#'
#' @description
#' Generates a random sample of y scores and x scores from a bivariate normal
#' distributions with specified population means, standard deviations, and 
#' correlation. This function is useful for generating hypothetical data for
#' classroom demonstrations.
#'
#'  
#' @param   n     sample size
#' @param   my    population mean of y scores
#' @param   mx    population mean of x scores
#' @param   sdy   population standard deviation of y scores
#' @param   sdx   population standard deviation of x scores
#' @param   cor   population correlation between x and y 
#' @param   dec   number of decimal points 
#'
#' 
#' @return 
#' Returns n pairs of y and x scores 
#' 
#' 
#' @examples
#' random.yx(10, 50, 20, 4, 2, .5, 1)
#'
#' # Should return: 
#' #        y    x
#' #  1  50.3 21.6
#' #  2  52.0 21.6
#' #  3  53.0 22.7
#' #  4  46.9 21.3
#' #  5  56.3 23.8
#' #  6  50.4 20.3
#' #  7  44.6 19.9
#' #  8  49.9 18.3
#' #  9  49.4 18.5
#' # 10  42.3 20.2
#'  
#' 
#' @export  
random.yx <- function(n, my, mx, sdy, sdx, cor, dec) {
 x0 <- rnorm(n, 0, 1)
 y0 <- cor*x0 + sqrt(1 - cor^2)*rnorm(n, 0, 1)
 x <- sdx*x0 + mx
 y <- sdy*y0 + my
 out <- as.data.frame(round(cbind(y, x), dec))
 colnames(out) <- c("y", "x")
 return(out)
}


#  sim.ci.cor ===============================================================
#' Simulates confidence interval coverage probability for a Pearson
#' correlation
#'
#'                               
#' @description
#' Performs a computer simulation of confidence interval performance for a 
#' Pearson correlation. A bias adjustment is used to reduce the bias of the
#' Fisher transformed Pearson correlation. Sample data can be generated 
#' from bivariate population distributions with five different marginal 
#' distributions. All distributions are scaled to have standard deviations
#' of 1.0. Bivariate random data with specified marginal skewness and
#' kurtosis are generated using the unonr function in the mnonr package. 
#'
#'
#' @param   alpha     alpha level for 1-alpha confidence
#' @param   n         sample size 
#' @param   cor       population Pearson correlation
#' @param   dist1     type of distribution for variable 1 (1, 2, 3, 4, or 5)
#' @param   dist2     type of distribution for variable 2 (1, 2, 3, 4, or 5)
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
#' * Coverage - probability of confidence interval including population correlation  
#' * Lower Error - probability of lower limit greater than population correlation
#' * Upper Error - probability of upper limit less than population correlation
#' * Ave CI Width - average confidence interval width
#'
#'
#' @examples
#' sim.ci.cor(.05, 30, .7, 4, 5, 1000)
#'
#' # Should return (within sampling error):
#' #      Coverage Lower Error Upper Error Ave CI Width
#' # [1,]  0.93815     0.05125      0.0106    0.7778518
#'
#'
#' @importFrom stats qt
#' @importFrom mnonr unonr
#' @export
sim.ci.cor <- function(alpha, n, cor, dist1, dist2, rep) {
 zcrit <- qnorm(1 - alpha/2)
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
#' Performs a computer simulation of confidence interval performance for a 
#' Spearman correlation. Sample data can be generated from bivariate population
#' distributions with five different marginal distributions. All distributions
#' are scaled to have standard deviations of 1.0. Bivariate random data with 
#' specified marginal skewness and kurtosis are generated using the unonr 
#' function in the mnonr package. 
#'
#'
#' @param   alpha     alpha level for 1-alpha confidence
#' @param   n         sample size 
#' @param   cor       population Spearman correlation
#' @param   dist1     type of distribution for variable 1 (1, 2, 3, 4, or 5)
#' @param   dist2     type of distribution for variable 2 (1, 2, 3, 4, or 5)
#' * 1 = Gaussian (skewness = 0 and excess kurtosis = 0) 
#' * 2 = platykurtic (skewness = 0 and excess kurtosis = -1.2)
#' * 3 = leptokurtic (skewness = 0 and excess kurtsois = 6)
#' * 4 = moderate skew (skewness = 1 and excess kurtosis = 1.5)
#' * 5 = large skew (skewness = 2 and excess kurtosis = 6)
#' @param   rep      number of Monte Carlo samples
#'  
#' 
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Coverage - probability of confidence interval including population correlation  
#' * Lower Error - probability of lower limit greater than population correlation
#' * Upper Error - probability of upper limit less than population corrrelation
#' * Ave CI Width - average confidence interval width
#'
#'
#' @examples
#' sim.ci.spear(.05, 30, .7, 4, 5, 1000)
#'
#' # Should return (within sampling error):
#' #      Coverage Lower Error Upper Error Ave CI Width
#' # [1,]  0.96235     0.01255      0.0251    0.4257299
#'
#'
#' @importFrom stats qt
#' @importFrom mnonr unonr
#' @export
sim.ci.spear <- function(alpha, n, cor, dist1, dist2, rep) {
 zcrit <- qnorm(1 - alpha/2)
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

