#  ======================= Confidence Intervals ==============================
#  ci.cor  ===================================================================
#' Confidence interval for a Pearson or partial correlation
#'
#'
#' @description
#' Computes a Fisher confidence interval for a population Pearson correlation  
#' or partial correlation with s control variables. Set s = 0 for a Pearson 
#' correlation. A bias adjustment is used to reduce the bias of the Fisher
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
#' * Estimate - estimated correlation (from input)
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' ci.cor(.05, .536, 0, 50)
#'
#' # Should return:
#' # Estimate        SE        LL        UL
#' #    0.536 0.1018149 0.2978573 0.7058914
#'  
#' 
#' @importFrom stats qnorm
#' @export
ci.cor <- function(alpha, cor, s, n) {
 if (cor > .999 || cor < -.999) {stop("correlation must be between -.999 and .999")}
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
 rownames(out) <- ""
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
#' * Estimate - estimated semipartial correlation (from input)
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
#' # Estimate        SE        LL        UL
#' #    0.582 0.1374298 0.2525662 0.7905182
#'  
#' 
#' @importFrom stats qnorm
#' @export
ci.spcor <- function(alpha, cor, r2, n) {
 if (cor > .999 || cor < -.999) {stop("correlation must be between -.999 and .999")}
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
 rownames(out) <- ""
 return(out)
}


#  ci.cor2 ====================================================================
#' Confidence interval for a 2-group Pearson correlation difference 
#'
#'
#' @description
#' Computes a confidence interval for a difference in population Pearson 
#' correlations in a 2-group design. A bias adjustment is used to 
#' reduce the bias of each Fisher transformed correlation.  
#'
#'  
#' @param  alpha	alpha level for 1-alpha confidence
#' @param  cor1	  	estimated Pearson correlation for group 1
#' @param  cor2	  	estimated Pearson correlation for group 2
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
#' # Estimate         SE         LL        UL
#' #    0.084 0.02967934 0.02803246 0.1463609
#'  
#' 
#' @importFrom stats qnorm
#' @export
ci.cor2 <- function(alpha, cor1, cor2, n1, n2) {
 if (cor1 > .999 || cor1 < -.999) {stop("correlation must be between -.999 and .999")}
 if (cor2 > .999 || cor2 < -.999) {stop("correlation must be between -.999 and .999")}
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
 rownames(out) <- ""
 return(out)
}


#  ci.cor.dep ================================================================
#' Confidence interval for a difference in dependent Pearson correlations
#'
#'
#' @description 
#' Computes a confidence interval for a difference in population Pearson 
#' correlations that are estimated from the same sample and have one 
#' variable in common. A bias adjustment is used to reduce the bias
#' of each Fisher transformed correlation. An approximate standard error
#' is recovered from the confidence interval.
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
#' ci.cor.dep(.05, .396, .179, .088, 166)
#'
#' # Should return:
#' # Estimate        SE         LL       UL
#' #    0.217 0.1026986 0.01323072 0.415802
#'  
#' 
#' @importFrom stats qnorm
#' @export 
ci.cor.dep <- function(alpha, cor1, cor2, cor12, n) {
 if (cor1 > .999 || cor1 < -.999) {stop("cor1 must be between -.999 and .999")}
 if (cor2 > .999 || cor2 < -.999) {stop("cor2 must be between -.999 and .999")}
 if (cor12 > .999 || cor12 < -.999) {stop("cor12 must be between -.999 and .999")}
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
 se <- (ul - ll)/(2*z)
 out <- t(c(diff, se, ll, ul))
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 rownames(out) <- ""
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
#' The correlations could also be correlations between two latent factors.
#' The function requires a point estimate and a 100(1 - alpha)% confidence
#' interval for each correlation as input. The confidence intervals for 
#' each correlation can be obtained using the ci.fisher function.
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
#' # Estimate   LL        UL
#' #      0.2 0.07 0.3220656
#'  
#' 
#' @importFrom stats qnorm
#' @export 
ci.cor2.gen <- function(cor1, ll1, ul1, cor2, ll2, ul2) {
 if (cor1 > .999 || cor1 < -.999) {stop("cor1 must be between -.999 and .999")}
 if (cor2 > .999 || cor2 < -.999) {stop("cor2 must be between -.999 and .999")}
 diff <- cor1 - cor2
 ll <- diff - sqrt((cor1 - ll1)^2 + (ul2 - cor2)^2)
 ul <- diff + sqrt((ul1 - cor1)^2 + (cor2 - ll2)^2)
 out <- t(c(diff, ll, ul))
 colnames(out) <- c("Estimate", "LL", "UL")
 rownames(out) <- ""
 return(out)
}


#  ci.pbcor ==================================================================
#' Confidence intervals for point-biserial correlations
#'
#'
#' @description 
#' Computes confidence intervals for two types of population point-biserial 
#' correlations. One type uses a weighted average of the group variances 
#' and is appropriate for nonexperimental designs with simple random sampling
#' (but not stratified random sampling). The other type uses an unweighted 
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
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' ci.pbcor(.05, 28.32, 21.48, 3.81, 3.09, 40, 40)
#'
#' # Should return:
#' #              Estimate         SE        LL        UL
#' # Weighted:   0.7065799 0.04890959 0.5885458 0.7854471
#' # Unweighted: 0.7020871 0.05018596 0.5808366 0.7828948
#'  
#' 
#' @importFrom stats qnorm
#' @export 
ci.pbcor <- function(alpha, m1, m2, sd1, sd2, n1, n2) {
 z <- qnorm(1 - alpha/2)
 df1 <- n1 - 1
 df2 <- n2 - 1
 n <- n1 + n2
 p <- n1/n
 s1 <- sqrt((df1*sd1^2 + df2*sd2^2)/(n - 2))
 s2 <- sqrt((sd1^2 + sd2^2)/2)
 d1 <- (m1 - m2)/s1
 d2 <- (m1 - m2)/s2
 k <- (n - 2)/(n*p*(1 - p))
 cor1 <- d1/sqrt(d1^2 + k)
 cor2 <- d2/sqrt(d2^2 + 4)
 se.d1 <- sqrt(d1^2*(1/n1 + 1/n2)/8 + (sd1^2/n1 + sd2^2/n2)/s1^2)
 se.d2 <- sqrt(d2^2*(sd1^4/df1 + sd2^4/df2)/(8*s2^4) + (sd1^2/df1 + sd2^2/df2)/s2^2) 
 se.cor1 <- (k/(d1^2 + k)^(3/2))*se.d1
 se.cor2 <- (4/(d2^2 + 4)^(3/2))*se.d2
 lld1 <- d1 - z*se.d1
 uld1 <- d1 + z*se.d1
 lld2 <- d2 - z*se.d2
 uld2 <- d2 + z*se.d2
 llc1 <- lld1/sqrt(lld1^2 + k)
 ulc1 <- uld1/sqrt(uld1^2 + k)
 llc2 <- lld2/sqrt(lld2^2 + 4)
 ulc2 <- uld2/sqrt(uld2^2 + 4)
 out1 <- t(c(cor1, se.cor1, llc1, ulc1))
 out2 <- t(c(cor2, se.cor2, llc2, ulc2))
 out <- rbind(out1, out2)
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 rownames(out) <- c("Weighted:", "Unweighted:")
 return(out)
}


#  ci.spear ==================================================================
#' Confidence interval for a Spearman correlation
#'
#'
#' @description
#' Computes a Fisher confidence interval for a population Spearman correlation. 
#' This function is not appropropriate for ordered categorical variables.
#'
#'  
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  y	  vector of y scores 
#' @param  x	  vector of x scores (paired with y)
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimated Spearman correlation
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
#' #  Estimate         SE        LL        UL
#' # 0.8699639 0.08241326 0.5840951 0.9638297
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
 rownames(out) <- ""
 return (out)
}
	
	
#  ci.spear2 ====================================================================
#' Confidence interval for a 2-group Spearman correlation difference 
#'
#'
#' @description
#' Computes a confidence interval for a difference of population Spearman 
#' correlations in a 2-group design. This function is not appropropriate 
#' for ordered categorical variables.
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
#' # Estimate         SE         LL        UL
#' #     0.06 0.08124926 -0.1003977 0.2185085     
#'  
#' 
#' @importFrom stats qnorm
#' @export
ci.spear2 <- function(alpha, cor1, cor2, n1, n2) {
 if (cor1 > .999 || cor1 < -.999) {stop("cor1 must be between -.999 and .999")}
 if (cor2 > .999 || cor2 < -.999) {stop("cor2 must be between -.999 and .999")}
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
 rownames(out) <- ""
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
#' * SE - standard error
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
#' # Estimate        SE       LL       UL
#' #   2.3744 0.3314752 1.751678 3.218499
#'  
#' 
#' @importFrom stats qnorm
#' @importFrom stats sd
#' @export
ci.mape <- function(alpha, res, s) {
 n <- length(res)
 df <- n - s - 1
 z <- qnorm(1 - alpha/2)
 c <- n/(n - (s + 2)/2)
 mape <- mean(abs(res))
 kur <- (sd(res)/mape)^2
 se <- sqrt((kur - 1)/df)
 ll <- exp(log(c*mape) - z*se)
 ul <- exp(log(c*mape) + z*se)
 se.mape <- c*mape*se
 out <- t(c(c*mape, se.mape, ll, ul))
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 rownames(out) <- ""
 return(out)
}


#  ci.ratio.mape2  ===================================================================
#' Confidence interval for a ratio of mean absolute prediction errors in a
#' 2-group design
#'
#'                      
#' @description
#' Computes a confidence interval for a ratio of population mean absolute 
#' prediction errors from a general linear model in two independent groups.
#' The number of predictor variables can differ across groups and the two
#' models can be non-nested. This function requires a vector of estimated 
#' residuals from each group. This function does not assume zero excess 
#' kurtosis but does assume symmetry in the population prediction errors for
#' the two models.
#'
#'  
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  res1   vector of residuals from group 1
#' @param  res2   vector of residuals from group 2
#' @param  s1	  number of predictor variables used in group 1
#' @param  s2	  number of predictor variables used in group 2
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * MAPE1 - bias adjusted mean absolute prediction error for group 1
#' * MAPE2 - bias adjusted mean absolute prediction error for group 2
#' * MAPE1/MAPE2 - ratio of bias adjusted mean absolute prediction errors
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' res1 <- c(-2.70, -2.69, -1.32, 1.02, 1.23, -1.46, 2.21, -2.10, 2.56, -3.02
#'         -1.55, 1.46, 4.02, 2.34)
#' res2 <- c(-0.71, -0.89, 0.72, -0.35, 0.33 -0.92, 2.37, 0.51, 0.68, -0.85,
#'         -0.15, 0.77, -1.52, 0.89, -0.29, -0.23, -0.94, 0.93, -0.31 -0.04)
#' ci.ratio.mape2(.05, res1, res2, 1, 1)
#'
#' # Should return:
#' #   MAPE1     MAPE2 MAPE1/MAPE2       LL       UL
#' # 2.58087 0.8327273    3.099298 1.917003 5.010761
#'  
#' 
#' @importFrom stats qnorm
#' @importFrom stats sd
#' @export
ci.ratio.mape2 <- function(alpha, res1, res2, s1, s2) {
 z <- qnorm(1 - alpha/2)
 n1 <- length(res1)
 df1 <- n1 - s1 - 1
 c1 <- n1/(n1 - (s1 + 2)/2)
 mape1 <- mean(abs(res1))
 kur1 <- (sd(res1)/mape1)^2
 se1 <- sqrt((kur1 - 1)/df1)
 n2 <- length(res2)
 df2 <- n2 - s2 - 1
 c2 <- n2/(n2 - (s2 + 2)/2)
 mape2 <- mean(abs(res2))
 kur2 <- (sd(res2)/mape2)^2
 se2 <- sqrt((kur2 - 1)/df2)
 se <- sqrt(se1^2 + se2^2)
 ll <- exp(log(c1*mape1) - log(c2*mape2) - z*se)
 ul <- exp(log(c1*mape1) - log(c2*mape2) + z*se)
 out <- t(c(c1*mape1, c2*mape2,(c1*mape1)/(c2*mape2), ll, ul))
 colnames(out) <- c("MAPE1", "MAPE2", "MAPE1/MAPE2", "LL", "UL")
 rownames(out) <- ""
 return(out)
}


#  ci.condslope  ==============================================================
#' Confidence intervals for conditional (simple) slopes in a linear model
#'
#'
#' @description
#' Computes confidence intervals and test statistics for population 
#' conditional slopes (simple slopes) in a general linear model that
#' includes a predictor variable (x1), a moderator variable (x2), and
#' a product predictor variable (x1*x2). Conditional slopes are computed  
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
#' * p - two-sided p-value
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
#' Computes a confidence interval and test statistic for a linear contrast
#' of population regression coefficients (e.g., a y-intercept or a slope
#' coefficient) across groups in a multiple group regression model. Equality
#' of error variances across groups is not assumed. A Satterthwaite adjustment
#' to the degrees of freedom is used to improve the accuracy of the confidence
#' interval. 
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
#' * p - two-sided p-value
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
#' # Estimate        SE        t      df          p        LL       UL
#' #    1.303 0.5085838 2.562016 78.8197 0.01231256 0.2906532 2.315347
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
 rownames(out) <- ""
 return(out)
}


#  ci.fisher ================================================================= 
#' Fisher confidence interval 
#'
#' 
#' @description
#' Computes a Fisher confidence interval for any type of correlation (e.g., 
#' Pearson, Spearman, Kendall-tau, tetrachoric, phi, partial, semipartial, 
#' etc.) or ordinal association such as gamma, Somers' d, or tau-b. The 
#' correlation could also be between two latent factors obtained
#' from a SEM analysis (the Fisher CI will be more accurate than the 
#' large-sample CI from a SEM analysis). The standard error can be a 
#' traditional standard error, a bootstrap standard error, or a robust 
#' standard error from a SEM analysis.
#'
#'
#' @param   alpha  alpha value for 1-alpha confidence
#' @param   cor    estimated correlation or association coefficient 
#' @param   se     standard error of correlation or association coefficient
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - correlation (from input) 
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#'
#'
#' @examples
#' ci.fisher(.05, .641, .052)
#'
#' # Should return:
#' # Estimate        LL        UL
#' #    0.641 0.5276396 0.7319293
#'
#'
#' @importFrom stats qnorm
#' @export
ci.fisher <- function(alpha, cor, se) {
 if (cor > .999 || cor < -.999) {stop("correlation must be between -.999 and .999")}
 z <- qnorm(1 - alpha/2)
 zr <- log((1 + cor)/(1 - cor))/2
 ll0 <- zr - z*se/(1 - cor^2)
 ul0 <- zr + z*se/(1 - cor^2)
 ll <- (exp(2*ll0) - 1)/(exp(2*ll0) + 1)
 ul <- (exp(2*ul0) - 1)/(exp(2*ul0) + 1)
 out <- t(c(cor, ll, ul))
 colnames(out) <- c("Estimate", "LL", "UL")
 return(out)
}


#  ci.indirect ===============================================================
#' Confidence interval for an indirect effect
#'
#'
#' @description
#' Computes a Monte Carlo confidence interval (500,000 trials) for a population
#' unstandardized or standardized indirect effect in a path model and a Sobel 
#' standard error. This function is not recommended for a standardized indirect 
#' if the standardized slopes are greater than .4 The Monte Carlo method is
#' general in that the slope estimates and standard errors do not need to be 
#' OLS estimates with homoscedastic standard errors. For example, LAD slope 
#' estimates and their standard errors, OLS slope estimates and 
#' heteroscedastic-consistent standard errors also could be used. In models 
#' with no direct effects, distribution-free Theil-Sen slope estimates with
#' recovered standard errors (see \link[statpsych]{ci.theil}) also could be used.
#'
#'  
#' @param  alpha  alpha level for 1-alpha confidence  
#' @param  b1     slope estimate for first path
#' @param  b2     slope estimate for second path
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
#' # Estimate       SE       LL       UL
#' #   4.7616 1.625282 2.178812 7.972262
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
 rownames(out) <- ""
 return(out)
}


#  ci.lc.glm ==================================================================
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
#' @param   q      vector of contrast coefficients
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimate of linear function 
#' * SE - standard error
#' * t - t test statistic 
#' * df - degrees of freedom
#' * p - two-sided p-value 
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
#' #
#' # Estimate        SE       t df           p       LL       UL
#' # 2.931111 0.4462518 6.56829  7 0.000313428 1.875893 3.986329
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
 rownames(out) <- ""
 return(out)
}


#  ci.lc.gen.bs ===============================================================
#' Confidence interval for a linear contrast of parameters in a between-subjects
#' design
#'
#'                                              
#' @description
#' Computes the estimate, standard error, and approximate confidence interval 
#' for a linear contrast of any type of parameter where each parameter value has
#' been estimated from a different sample. The parameter values are assumed to 
#' be of the same type and their sampling distributions are assumed to be 
#' approximately normal.
#'
#'
#' @param  alpha   alpha level for 1-alpha confidence
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
#' # Estimate        SE       LL       UL
#' #     1.63 0.2573806 1.125543 2.134457
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
 rownames(out) <- ""
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
#' * R-squared - estimate of unadjusted R-squared (from input)
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
#' # R-squared adj R-squared         SE         LL        UL  
#' #     0.241     0.2206696 0.06752263 0.09819599 0.3628798
#'  
#' 
#' @importFrom stats qf
#' @export
ci.rsqr <- function(alpha, r2, s, n) {
 if (r2 > .999 || r2 < .001) {stop("squared multiple correlation must be between .001 and .999")}
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
 if (v1 < 1) {v1 = 1}
 if (v2 < 1) {v2 = 1}
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
   if (v1 < 1) {v1 = 1}
   if (v2 < 1) {v2 = 1}
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
 rownames(out) <- ""
 return(out)
}


#  ci.theil ===================================================================
#' Theil-Sen estimate and confidence interval for slope
#'
#'                                 
#' @description
#' Computes a Theil-Sen estimate and distribution-free confidence interval 
#' for the population slope in a simple linear regression model. An approximate 
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
#' # Estimate        SE        LL   UL
#' #      0.5 0.1085927 0.3243243 0.75
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
  rownames(out) <- ""
  return(out)
}


#  ci.cronbach2 ===============================================================
#' Confidence interval for a difference in Cronbach reliabilities in a 2-group 
#' design
#'
#'
#' @description
#' Computes a confidence interval for a difference in population Cronbach 
#' reliability coefficients in a 2-group design. The number of measurements
#' (e.g., items or raters) used in each group need not be equal.
#'
#'  
#' @param  alpha    alpha level for 1-alpha confidence
#' @param  rel1     estimated Cronbach reliability for group 1
#' @param  rel2     estimated Cronbach reliability for group 2
#' @param  r1       number of measurements used in group 1
#' @param  r2       number of measurements used in group 2
#' @param  n1       sample size for group 1
#' @param  n2       sample size for group 2
#'
#' 
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimated reliability difference
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @references
#' \insertRef{Bonett2015}{statpsych}
#'
#'
#' @examples
#' ci.cronbach2(.05, .88, .76, 8, 8, 200, 250)
#'
#' # Should return:
#' # Estimate         LL       UL
#' #     0.12 0.06973411 0.173236
#'  
#' 
#' @importFrom stats qf
#' @export
ci.cronbach2 <- function(alpha, rel1, rel2, r1, r2, n1, n2) {
 if (rel1 > .999 || rel1 < .001) {stop("rel1 must be between .001 and .999")}
 if (rel2 > .999 || rel2 < .001) {stop("rel2 must be between .001 and .999")}
 df11 <- n1 - 1
 df21 <- n1*(r1 - 1)
 f11 <- qf(1 - alpha/2, df11, df21)
 f21 <- qf(1 - alpha/2, df21, df11)
 f01 <- 1/(1 - rel1)
 LL1 <- 1 - f11/f01
 UL1 <- 1 - 1/(f01*f21)
 df12 <- n2 - 1
 df22 <- n2*(r2 - 1)
 f12 <- qf(1 - alpha/2, df12, df22)
 f22 <- qf(1 - alpha/2, df22, df12)
 f02 <- 1/(1 - rel2)
 LL2 <- 1 - f12/f02
 UL2 <- 1 - 1/(f02*f22)
 d <- rel1 - rel2
 ll <- d - sqrt((rel1 - LL1)^2 + (UL2 - rel2)^2)
 ul <- d + sqrt((UL1 - rel1)^2 + (rel2 - LL2)^2)
 out <- t(c(d, ll, ul))
 colnames(out) <- c("Estimate", "LL", "UL")
 rownames(out) <- ""
 return(out)
}


#  ci.rel2 ================================================================
#' Confidence interval for a 2-group reliability difference
#'
#'
#' @description 
#' Computes a 100(1 - alpha)% confidence interval for a difference in 
#' population reliabilities in a 2-group design. This function can be
#' used with any type of reliability coefficient (e.g., Cronbach alpha,
#' McDonald omega, intraclass reliability). The function requires a
#' point estimate and a 100(1 - alpha)% confidence interval for each
#' reliability as input. 
#'
#'  
#' @param  rel1  estimated reliability for group 1 
#' @param  ll1   lower limit for group 1 reliability
#' @param  ul1   upper limit for group 1 reliability
#' @param  rel2  estimated reliability for group 2
#' @param  ll2   lower limit for group 2 reliability
#' @param  ul2   upper limit for group 2 reliability
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimated reliability difference
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @references
#' \insertRef{Bonett2015}{statpsych}
#'
#'
#' @examples
#' ci.rel2(.4, .35, .47, .2, .1, .32)
#'
#' # Should return:
#' # Estimate   LL        UL
#' #      0.2 0.07 0.3220656
#'  
#' 
#' @importFrom stats qnorm
#' @export 
ci.rel2 <- function(rel1, ll1, ul1, rel2, ll2, ul2) {
 if (rel1 > .999 || rel1 < .001) {stop("reliability must be between .001 and .999")}
 if (rel2 > .999 || rel2 < .001) {stop("reliability must be between .001 and .999")}
 diff <- rel1 - rel2
 ll <- diff - sqrt((rel1 - ll1)^2 + (ul2 - rel2)^2)
 ul <- diff + sqrt((ul1 - rel1)^2 + (rel2 - ll2)^2)
 out <- t(c(diff, ll, ul))
 colnames(out) <- c("Estimate", "LL", "UL")
 rownames(out) <- ""
 return(out)
}


#  ci.bscor ==================================================================
#' Confidence interval for a biserial correlation
#'
#'
#' @description 
#' Computes a confidence interval for a population biserial correlation. A
#' biserial correlation can be used when one variable is quantitative and the 
#' other variable has been artificially dichotomized to create two groups.
#' The biserial correlation estimates the correlation between the observed 
#' quantitative variable and the unobserved quantitative variable that has 
#' been measured on a dichotomous scale. 
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
#' @details
#' This function computes a point-biserial correlation and its standard error
#' as a function of a standardized mean difference with a weighted variance
#' standardizer. Then the point-biserial estimate is transformed into a 
#' biserial correlation using the traditional adjustment. The adjustment is 
#' also applied to the point-biserial standard error to obtain the standard 
#' error for the biserial correlation. 
#' 
#' The biserial correlation assumes that the observed quantitative variable 
#' and the unobserved quantitative variable have a bivariate normal 
#' distribution. Bivariate normality is a crucial assumption underlying the
#' transformation of a point-biserial correlation to a biserial correlation.
#' Bivariate normality also implies equal variances of the observed 
#' quantitative variable at each level of the dichotomized variable, and this
#' assumption is made in the computation of the standard error.
#'
#'
#' @references
#' \insertRef{Bonett2020a}{statpsych}
#'
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimated biserial correlation
#' * SE - standard error
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' ci.bscor(.05, 28.32, 21.48, 3.81, 3.09, 40, 40)
#'
#' # Should return:
#' #   Estimate         SE        LL        UL
#' #  0.8855666 0.06129908  0.7376327 0.984412
#'  
#' 
#' @importFrom stats qnorm
#' @importFrom stats dnorm
#' @export 
ci.bscor <- function(alpha, m1, m2, sd1, sd2, n1, n2) {
 z <- qnorm(1 - alpha/2)
 df1 <- n1 - 1
 df2 <- n2 - 1
 u <- n1/(n1 + n2)
 a <- sqrt(u*(1 - u))/dnorm(qnorm(u))
 s <- sqrt((df1*sd1^2 + df2*sd2^2)/(df1 + df2))
 d <- (m1 - m2)/s
 c <- (df1 + df2)/((n1 + n2)*(u*(1 - u)))
 pbcor <- d/sqrt(d^2 + c)
 bscor <- pbcor*a
 if (bscor > 1) {bscor = .99999}
 if (bscor < -1) {bscor = -.99999}
 se.d <- sqrt(d^2*(1/n1 + 1/n2)/8 + 1/n1 + 1/n2)
 se.pbcor <- (c/(d^2 + c)^(3/2))*se.d  
 se.bscor <- se.pbcor*a 
 lld <- d - z*se.d
 uld <- d + z*se.d
 ll <- a*lld/sqrt(lld^2 + c)
 ul <- a*uld/sqrt(uld^2 + c)
 if (ul > 1) {ul = 1}
 if (ll < -1) {ll = -1}
 out <- t(c(bscor, se.bscor, ll, ul))
 colnames(out) <- c("Estimate", "SE", "LL", "UL")
 rownames(out) <- ""
 return(out)
}


#  pi.cor ===================================================================== 
#' Prediction limits for an estimated correlation
#'
#'                                        
#' @description
#' Computes approximate one-sided or two-sided prediction limits for the 
#' estimated Pearson correlation in a future study with a planned sample 
#' size of n. The prediction interval uses a correlation estimate from a
#' prior study that had a sample size of n0. 
#'
#' Several confidence interval sample size functions in this package require
#' a planning value of the estimated Pearson correlation that is expected 
#' in the planned study. A one-sided lower correlation prediction limit is 
#' useful as a correlation planning value for the sample size required to 
#' obtain a confidence interval with desired width. This strategy for 
#' specifying a correlation planning value is useful in applications where 
#' the population correlation in the prior study is assumed to be very similar
#' to the population correlation in the planned study. 
#'
#'
#' @param  alpha  alpha value for 1-alpha confidence 
#' @param  cor    estimated Pearson correlation from prior study
#' @param  n0     sample size used to estimate correlation in prior study
#' @param  n      planned sample size of future study
#' @param  type   
#' * set to 1 for two-sided prediction interval 
#' * set to 2 for one-sided upper prediction limit 
#' * set to 3 for one-sided lower prediction limit 
#'
#'
#' @return 
#' Returns one-sided or two-sided prediction limit(s) of an estimated 
#' Pearson correlation in a future study
#'
#'
#' @examples
#' pi.cor(.1, .761, 50, 100, 1)
#'
#' # Should return:
#' #         LL        UL
#' #  0.6034092 0.8573224
#'  
#' pi.cor(.1, .761, 50, 100, 3)
#'
#' # Should return:
#' #         LL
#' #  0.6428751
#'  
#' 
#' @importFrom stats qnorm
#' @export
pi.cor <- function(alpha, cor, n0, n, type) {
 if (type == 1) {
  z <- qnorm(1 - alpha/2)
  cor.z <- log((1 + cor)/(1 - cor))/2
  ll0 <- cor.z - cor/(2*(n0 - 1)) - z*sqrt(1/(n0 - 3) + 1/(n - 3))
  ul0 <- cor.z - cor/(2*(n0 - 1)) + z*sqrt(1/(n0 - 3) + 1/(n - 3))
  ll <- (exp(2*ll0) - 1)/(exp(2*ll0) + 1)
  ul <- (exp(2*ul0) - 1)/(exp(2*ul0) + 1)
  out <- t(c(ll, ul))
  colnames(out) <- c("LL", "UL")
 }
 else if (type == 2) {
  z <- qnorm(1 - alpha)
  cor.z <- log((1 + cor)/(1 - cor))/2
  ul0 <- cor.z - cor/(2*(n0 - 1)) + z*sqrt(1/(n0 - 3) + 1/(n - 3))
  ul <- (exp(2*ul0) - 1)/(exp(2*ul0) + 1)
  out <- matrix(ul, nrow = 1, ncol = 1)
  colnames(out) <- "UL"
 }
 else {
  z <- qnorm(1 - alpha)
  cor.z <- log((1 + cor)/(1 - cor))/2
  ll0 <- cor.z - cor/(2*(n0 - 1)) - z*sqrt(1/(n0 - 3) + 1/(n - 3))
  ll <- (exp(2*ll0) - 1)/(exp(2*ll0) + 1)
  out <- matrix(ll, nrow = 1, ncol = 1)
  colnames(out) <- "LL"
 }
 rownames(out) <- ""
 return(out)
}


# ci.bayes.cor ============================================================
#' Bayesian credible interval for a Pearson or partial correlation with a
#' skeptical prior
#'
#'
#' @description
#' Computes an approximate Bayesian credible interval for a Pearson or  
#' partial correlation with a skeptical prior. The skeptical prior 
#' distribution is Normal with a mean of 0 and a small standard deviation.
#' A skeptical prior assumes that the population correlation is within 
#' a range of small values (-r to r). If the skeptic is 95% confident that
#' the population correlation is between -r and r, then the prior standard
#' deviation can be set to r/1.96. A correlation that is less than .2 in 
#' absolute value is typically considered to be "small", and the prior 
#' standard deviation could then be set to .2/1.96. A correlation value
#' that is considered to be small will depend on the application. Set s = 0
#' for a Pearson correlation. 
#'
#'
#' @param   alpha        alpha level for 1-alpha credibility interval
#' @param   prior_sd     standard deviation of skeptical prior distribution 
#' @param   cor          estimated Pearson or partial correlation
#' @param   s	     	 number of control variables
#' @param   n            sample size
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Posterior mean - posterior mean (Bayesian estimate of correlation)
#' * LL - lower limit of the credible interval
#' * UL - upper limit of the credible interval
#'
#'
#' @examples
#' ci.bayes.cor(.05, .1, .536, 0, 50)
#'
#' # Should return:
#' # Posterior mean         LL        UL
#' #      0.1873765 0.02795441 0.3375031
#'
#' ci.bayes.cor(.05, .1, .536, 0, 300)
#'
#' # Should return:
#' #  Posterior mean        LL        UL
#' #       0.4195068 0.3352449 0.4971107
#'
#'
#' @importFrom stats qnorm
#' @export
ci.bayes.cor <- function(alpha, prior_sd, cor, s, n) {
 z <- qnorm(1 - alpha/2)
 se <- 1/sqrt(n - s - 3)
 zr <- log((1 + cor)/(1 - cor))/2 - cor/(2*(n - 1))
 post_sd <- sqrt(1/(1/prior_sd^2 + 1/se^2))
 post_mean <- (zr/se^2)/(1/prior_sd^2 + 1/se^2)
 ll0 <- post_mean - z*post_sd
 ul0 <- post_mean + z*post_sd
 mean <- (exp(2*post_mean) - 1)/(exp(2*post_mean) + 1)
 ll <- (exp(2*ll0) - 1)/(exp(2*ll0) + 1)
 ul <- (exp(2*ul0) - 1)/(exp(2*ul0) + 1)
 out <- t(c(mean, ll, ul))
 colnames(out) <- c("Posterior mean", "LL", "UL")
 rownames(out) <- ""
 return(out)
}


# ci.bayes.spcor ============================================================
#' Bayesian credible interval for a semipartial correlation with a
#' skeptical prior
#'
#'
#' @description
#' Computes an approximate Bayesian credible interval for a semipartial 
#' correlation with a skeptical prior. The skeptical prior distribution is
#' Normal with a mean of 0 and a small standard deviation. A skeptical prior 
#' assumes that the population semipartial correlation is within a range of 
#' small values (-r to r). If the skeptic is 95% confident that the population
#' correlation is between -r and r, then the prior standard deviation can be 
#' set to r/1.96. A semipartial correlation that is less than .2 in absolute 
#' value is typically considered to be "small", and the prior standard 
#' deviation could then be set to .2/1.96. A semipartial correlation value
#' that is considered to be small will depend on the application. This function 
#' requires the standard error of the estimated semipartial correlation which
#' can be obtained from the ci.spcor function. 
#'
#'
#' @param   alpha        alpha level for 1-alpha credibility interval
#' @param   prior_sd     standard deviation of skeptical prior distribution 
#' @param   cor          estimated semipartial partial correlation
#' @param   se	     	 standard error of estimated semipartial correlation
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Posterior mean - posterior mean (Bayesian estimate of correlation)
#' * LL - lower limit of the credible interval
#' * UL - upper limit of the credible interval
#'
#'
#' @examples
#' ci.bayes.spcor(.05, .1, .582, .137)
#'
#' # Should return:
#' #  Posterior mean        LL        UL
#' #       0.2272797 0.07288039 0.3710398
#'
#'
#' @importFrom stats qnorm
#' @export
ci.bayes.spcor <- function(alpha, prior_sd, cor, se) {
 z <- qnorm(1 - alpha/2)
 zr <- log((1 + cor)/(1 - cor))/2 
 post_sd <- sqrt(1/(1/prior_sd^2 + 1/se^2))
 post_mean <- (zr/se^2)/(1/prior_sd^2 + 1/se^2)
 ll0 <- post_mean - z*post_sd
 ul0 <- post_mean + z*post_sd
 mean <- (exp(2*post_mean) - 1)/(exp(2*post_mean) + 1)
 ll <- (exp(2*ll0) - 1)/(exp(2*ll0) + 1)
 ul <- (exp(2*ul0) - 1)/(exp(2*ul0) + 1)
 out <- t(c(mean, ll, ul))
 colnames(out) <- c("Posterior mean", "LL", "UL")
 rownames(out) <- ""
 return(out)
}


#  ci.slope.mean.bs ===========================================================
#' Confidence interval for the slope of means in a one-factor experimental 
#' design with a quantitative between-subjects factor
#' 
#' 
#' @description
#' Computes a test statistic and confidence interval for the slope of means in 
#' a one-factor experimental design with a quantitative between-subjects 
#' factor. This function computes both the unequal variance and equal variance
#' confidence intervals and test statistics. A Satterthwaite adjustment to the
#' degrees of freedom is used with the unequal variance method. 
#'
#'
#' @param     alpha  	alpha level for 1-alpha confidence
#' @param     m     	vector of sample means
#' @param     sd    	vector of sample standard deviations
#' @param     n     	vector of sample sizes
#' @param     x     	vector of quantiative factor values
#' 
#'
#' @return 
#' Returns a 2-row matrix. The columns are:
#' * Estimate - estimated slope
#' * SE - standard error
#' * t - t test statistic
#' * df - degrees of freedom
#' * p - two-sided p-value
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


#  ci.slope.median.bs =========================================================
#' Confidence interval for the slope of medians in a one-factor experimental 
#' design with a quantitative between-subjects factor
#' 
#' 
#' @description
#' Computes a distrbution-free test and confidence interval for the slope 
#' of medians in a one-factor experimental design with a quantitative 
#' between-subjects factor using sample group medians and standard errors
#' as input. The sample median and standard error for each group can be 
#' computed using the \link[statpsych]{ci.median} function. 
#'
#'
#' @param     alpha  	alpha level for 1-alpha confidence
#' @param     m     	vector of sample median
#' @param     se    	vector of standard errors
#' @param     x     	vector of quantitative factor values
#' 
#'
#' @return 
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimated slope
#' * SE - standard error
#' * z - z test statistic 
#' * p - two-sided p-value
#' * LL - lower limit of the confidence interval
#' * UL - upper limit of the confidence interval
#' 
#' 
#' @examples
#' m <- c(33.5, 37.9, 38.0, 44.1)
#' se <- c(0.84, 0.94, 1.65, 2.98)
#' x <- c(5, 10, 20, 30)
#' ci.slope.median.bs(.05, m, se, x)
#'
#' # Should return:
#' #   Estimate        SE        z           p        LL        UL
#' #  0.3664407 0.1163593 3.149216 0.001637091 0.1383806 0.5945008
#'
#'
#' @importFrom stats qnorm
#' @importFrom stats pnorm
#' @export
ci.slope.median.bs <- function(alpha, m, se, x) {
 zcrit <- qnorm( 1 - alpha/2)
 mx <- mean(x)
 ssx <- sum((x - mx)^2)
 v <- (x - mx)/ssx
 est <- t(v)%*%m 
 se <- sqrt(t(v)%*%diag(se^2)%*%v)
 z <- est/se
 p <- 2*(1 - pnorm(abs(z)))
 ll <- est - zcrit*se
 ul <- est + zcrit*se
 out <- t(c(est, se, z, p, ll, ul))
 colnames(out) <- c("Estimate", "SE", "z", "p", "LL", "UL")
 rownames(out) <- ""
 return(out)
}

#  ======================== Hypothesis Tests ==================================
# test.cor ===================================================================
#' Hypothesis test for a Pearson or partial correlation 
#'
#'                        
#' @description
#' Computes a t test for a test of the null hypothesis that a population 
#' Pearson or partial correlations is equal to 0, or a z test using a Fisher 
#' transformation for a test of the null hypothesis that a Pearson or
#' partial correlation is equal to some specified nonzero value. Set s = 0 
#' for a Pearson correlation. The hypothesis testing results should be 
#' accompanied with a confidence interval for the population Pearson or
#' partial correlation value (see \link[statpsych]{ci.cor}).
#'
#'
#' @param  cor     estimated correlation 
#' @param  n       sample size 
#' @param  s       number of control variables
#' @param  h       null hypothesis value of correlation
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimate of correlation 
#' * t or z - t test statistic (for h = 0) or z test statistic (for nonzero h)
#' * p - two-sided p-value
#'
#'
#' @examples
#' test.cor(.484, 100, 0, .2)
#'
#' # Should return:
#' # Estimate        z           p
#' #    0.484 3.205432 0.001348601
#'
#'
#' test.cor(.372, 100, 0, 0)
#'
#' # Should return:
#' #  Estimate        t df           p
#' #     0.372 3.967337 98 0.000138436
#'
#'
#' @importFrom stats pnorm
#' @export
test.cor <- function(cor, n, s, h) {
 if (cor > .9999) {stop("correlation cannot be greater than .9999")}
 if (cor < -.9999) {stop("correlation cannot be less than -.9999")}
 if (h == 0) {
  df <- n - 2 - s
  se <- sqrt((1 - cor^2)/df)
  t <- cor/se
  pval <- 2*(1 - pt(abs(t), df))
  out <- t(c(cor, t, df, pval))
  colnames(out) <- c("Estimate", "t", "df", "p")
 }
 else {
   r.z <- log((1 + cor)/(1 - cor))/2
   h.z <- log((1 + h)/(1 - h))/2
   se.z <- sqrt(1/(n - 3 - s))
   z <- (r.z - h.z)/se.z
   pval <- 2*(1 - pnorm(abs(z)))
   out <- t(c(cor, z, pval))
   colnames(out) <- c("Estimate", "z", "p")
 }
 rownames(out) <- ""
 return(out)
}


# test.spear ===================================================================
#' Hypothesis test for a Spearman correlation 
#'
#'                      
#' @description
#' Computes a t test for a test of the null hypothesis that a population 
#' Spearman correlation is equal to 0, or a z test using a Fisher transformation
#' for a test of the null hypothesis that a Spearman correlation is equal to
#' some specified nonzero value. The hypothesis testing results should be 
#' accompanied with a confidence interval for the population Spearman
#' correlation value (see \link[statpsych]{ci.spear}).
#'
#'
#' @param  cor     estimated correlation 
#' @param  h       null hypothesis value of correlation
#' @param  n       sample size 
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimate of correlation 
#' * t or z - t test statistic (for h = 0) or z test statistic (for nonzero h)
#' * p - two-sided p-value
#'
#'
#' @examples
#' test.spear(.471, .2, 100)
#'
#' # Should return:
#' # Estimate        z           p
#' #     0.471 3.009628 0.00261568
#'
#'
#' test.spear(.342, 0, 100)
#'
#' # Should return:
#' #  Estimate        t df            p
#' #     0.342 3.602881 98 0.0004965008
#'
#'
#' @importFrom stats pnorm
#' @export
test.spear <- function(cor, h, n) {
 if (cor > .9999) {stop("correlation cannot be greater than .9999")}
 if (cor < -.9999) {stop("correlation cannot be less than -.9999")}
 if (h == 0) {
  df <- n - 2
  se <- sqrt((1 - cor^2)/df)
  t <- cor/se
  pval <- 2*(1 - pt(abs(t), df))
  out <- t(c(cor, t, df, pval))
  colnames(out) <- c("Estimate", "t", "df", "p")
 }
 else {
   r.z <- log((1 + cor)/(1 - cor))/2
   h.z <- log((1 + h)/(1 - h))/2
   se.z <- sqrt((1 + h^2/2)/(n - 3))
   z <- (r.z - h.z)/se.z
   pval <- 2*(1 - pnorm(abs(z)))
   out <- t(c(cor, z, pval))
   colnames(out) <- c("Estimate", "z", "p")
 }
 rownames(out) <- ""
 return(out)
}


# test.cor2 =================================================================
#' Hypothesis test for a 2-group Pearson or partial correlation difference
#'
#'
#' @description
#' Computes a z test for a difference of population Pearson or partial 
#' correlations in a 2-group design. Set s = 0 for a Pearson correlation. 
#' The hypothesis testing results should be accompanied with a confidence 
#' interval for the difference in population correlation values 
#' (see \link[statpsych]{ci.cor2}).
#'
#'
#' @param  cor1    estimated correlation for group 1
#' @param  cor2    estimated correlation for group 2
#' @param  n1      sample size for group 1
#' @param  n2      sample size for group 2
#' @param  s       number of control variables
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimate of correlation difference
#' * z - z test statistic
#' * p - two-sided p-value
#'
#'
#' @examples
#' test.cor2(.684, .437, 100, 125, 0)
#'
#' # Should return:
#' # Estimate        z           p
#' #    0.247 2.705709 0.006815877
#'
#'
#' @importFrom stats pnorm
#' @export
test.cor2 <- function(cor1, cor2, n1, n2, s) {
 if (cor1 > .9999) {stop("correlation cannot be greater than .9999")}
 if (cor1 < -.9999) {stop("correlation cannot be less than -.9999")}
 if (cor2 > .9999) {stop("correlation cannot be greater than .9999")}
 if (cor2 < -.9999) {stop("correlation cannot be less than -.9999")}
 z1 <- log((1 + cor1)/(1 - cor1))/2
 z2 <- log((1 + cor2)/(1 - cor2))/2
 v1 <- 1/(n1 - 3 - s)
 v2 <- 1/(n2 - 3 - s)
 se <- sqrt(v1 + v2)
 diff <- cor1 - cor2
 z <- (z1 - z2)/se
 pval <- 2*(1 - pnorm(abs(z)))
 out <- t(c(diff, z, pval))
 colnames(out) <- c("Estimate", "z", "p")
 rownames(out) <- ""
 return(out)
}


# test.spear2 =================================================================
#' Hypothesis test for a 2-group Spearman correlation difference
#'
#'
#' @description
#' Computes a z test for a difference of population Spearman correlations in a 
#' 2-group design. The test statistic uses a Bonett-Wright standard error for 
#' each Spearman correlation. The hypothesis testing results should be 
#' accompanied with a confidence interval for a difference in population 
#' Spearman correlation values (see \link[statpsych]{ci.spear2}).
#'
#'
#' @param  cor1    estimated Spearman correlation for group 1
#' @param  cor2    estimated Spearman correlation for group 2
#' @param  n1      sample size for group 1
#' @param  n2      sample size for group 2
#'
#'
#' @return
#' Returns a 1-row matrix. The columns are:
#' * Estimate - estimate of correlation difference
#' * z - z test statistic
#' * p - two-sided p-value
#'
#'
#' @references
#' \insertRef{Bonett2000}{statpsych}
#'
#'
#' @examples
#' test.spear2(.684, .437, 100, 125)
#'
#' # Should return:
#' # Estimate        z          p
#' #    0.247 2.498645 0.01246691
#'
#'
#' @importFrom stats pnorm
#' @export
test.spear2 <- function(cor1, cor2, n1, n2) {
 if (cor1 > .9999) {stop("correlation cannot be greater than .9999")}
 if (cor1 < -.9999) {stop("correlation cannot be less than -.9999")}
 if (cor2 > .9999) {stop("correlation cannot be greater than .9999")}
 if (cor2 < -.9999) {stop("correlation cannot be less than -.9999")}
 z1 <- log((1 + cor1)/(1 - cor1))/2
 z2 <- log((1 + cor2)/(1 - cor2))/2
 v1 <- (1 + cor1^2/2)/(n1 - 3)
 v2 <- (1 + cor2^2/2)/(n2 - 3)
 se <- sqrt(v1 + v2)
 diff <- cor1 - cor2
 z <- (z1 - z2)/se
 pval <- 2*(1 - pnorm(abs(z)))
 out <- t(c(diff, z, pval))
 colnames(out) <- c("Estimate", "z", "p")
 rownames(out) <- ""
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
#' that the population means are monotonic decreasing. If one or more upper 
#' limits are less than 0 and no lower limits are greater than 0, then
#' conclude that the population means are monotonic increasing. Reject the 
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


#  test.mono.median.bs ==========================================================
#' Test of a monotonic trend in medians for an ordered between-subjects factor
#' 
#'                     
#' @description
#' Computes simultaneous confidence intervals for all adjacent pairwise
#' comparisons of population medians using sample group medians and 
#' standard errors as input. If one or more lower limits are greater than 0  
#' and no upper limit is less than 0, then conclude that the population 
#' medians are monotonic decreasing. If one or more upper limits are less 
#' than 0 and no lower limits are greater than 0, then conclude that the 
#' population medians are monotonic increasing. Reject the hypothesis of a
#' monotonic trend if any lower limit is greater than 0 and any upper limit
#' is less than 0. The sample median and standard error for each group
#' can be computed using the \link[statpsych]{ci.median} function. 
#'
#'
#' @param  alpha   alpha level for simultaneous 1-alpha confidence
#' @param  m       vector of estimated group medians
#' @param  se      vector of estimated group standard errors
#'
#'
#' @return 
#' Returns a matrix with the number of rows equal to the number
#' of adjacent pairwise comparisons. The columns are:
#' * Estimate - estimated median difference
#' * SE - standard error
#' * LL - one-sided lower limit of the confidence interval
#' * UL - one-sided upper limit of the confidence interval
#'
#'
#' @examples
#' m <- c(12.86, 24.57, 36.29, 53.21)
#' se <- c(2.85, 2.99, 3.73, 3.88)
#' test.mono.median.bs(.05, m, se)
#'
#' # Should return:
#' #      Estimate       SE        LL         UL
#' #  1 2   -11.71 4.130690 -21.59879 -1.8212115
#' #  2 3   -11.72 4.780481 -23.16438 -0.2756247
#' #  3 4   -16.92 5.382128 -29.80471 -4.0352947
#'
#'
#' @importFrom stats qnorm
#' @export
test.mono.median.bs <-function(alpha, m, se) {
 a <- length(m)
 v <- se^2
 m1 <- m[1: a - 1]
 m2 <- m[2: a]
 Estimate <- m1 - m2
 v1 <- v[1: a - 1]
 v2 <- v[2: a]
 SE <- sqrt(v1 + v2)
 zcrit <- qnorm(1 - alpha/(2*(a - 1)))
 LL <- Estimate - zcrit*SE
 UL <- Estimate + zcrit*SE
 pair = cbind(seq(1, a - 1), seq(2, a))
 out <- cbind(pair, Estimate, SE, LL, UL)
 rownames(out) <- rep("", a - 1)
 return(out)
}


#  =================== Sample Size for Desired Precision ======================
#  size.ci.slope ==============================================================
#' Sample size for a slope confidence interval
#'
#'
#' @description
#' Computes the total sample size required to estimate a population slope with
#' desired confidence interval precision in a between-subjects design with a 
#' quantitative factor. In an experimental design, the total sample size 
#' would be allocated to the levels of the quantitative factor and it might
#' be necessary to increase the total sample size to achieve equal sample
#' sizes. Set the error variance planning value to the largest value within 
#' a plausible range for a conservatively large sample size.
#'
#'  
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  evar   planning value of within-group (error) variance
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
#' # Total sample size
#' #                83
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
 colnames(out) <- "Total sample size"
 rownames(out) <- ""
 return(out)
}


#  size.ci.cor ==============================================================
#' Sample size for a Pearson or partial correlation confidence interval 
#'
#'
#' @description
#' Computes the sample size required to estimate a population Pearson or
#' partial correlation with desired confidence interval precision. 
#' Set s = 0 for a Pearson correlation. Set the correlation planning value
#' to the smallest absolute value within a plausible range for a conservatively 
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
#' # Sample size
#' #         188
#'  
#' 
#' @importFrom stats qnorm
#' @export  
size.ci.cor <- function(alpha, cor, s, w) {
 if (cor > .999 | cor < -.999) {stop("correlation must be between -.999 and .999")}
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
 rownames(out) <- ""
 return(out)
}


#  size.ci.spear ==============================================================
#' Sample size for a Spearman correlation confidence interval 
#'
#'
#' @description
#' Computes the sample size required to estimate a population Spearman correlation
#' with desired confidence interval precision. Set the correlation planning value
#' to the smallest absolute value within a plausible range for a conservatively
#' large sample size.
#'
#'  
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  cor    planning value of Spearman correlation
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
#' size.ci.spear(.05, .362, .25)
#'
#' # Should return:
#' # Sample size
#' #         200
#'  
#' 
#' @importFrom stats qnorm
#' @export  
size.ci.spear <- function(alpha, cor, w) {
 if (cor > .999 | cor < -.999) {stop("correlation must be between -.999 and .999")}
 z <- qnorm(1 - alpha/2)
 n1 <- ceiling(4*(1 + cor^2/2)*(1 - cor^2)^2*(z/w)^2 + 3)
 zr <- log((1 + cor)/(1 - cor))/2 
 se <- sqrt((1 + cor^2/2)/(n1 - 3))
 ll0 <- zr - z*se
 ul0 <- zr + z*se
 ll <- (exp(2*ll0) - 1)/(exp(2*ll0) + 1)
 ul <- (exp(2*ul0) - 1)/(exp(2*ul0) + 1)
 n <- ceiling((n1 - 3)*((ul - ll)/w)^2 + 3)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 rownames(out) <- ""
 return(out)
}


#  size.ci.pbcor ==============================================================
#' Sample size for a point-biserial correlation confidence interval 
#'
#'
#' @description
#' Computes the sample size required to estimate a population point-biserial 
#' correlation with desired confidence interval precision in a two-group 
#' nonexperimental design with simple random sampling. A two-group 
#' nonexperimental design implies two subpopulations (e.g., all boys and all
#' girls in a school district). This function requires a planning value for the
#' proportion of population members who belong to one of the two subpopulations. 
#' Set the correlation planning value to the smallest absolute value within a 
#' plausible range for a conservatively large sample size.
#'
#'  
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  cor    planning value of point-biserial correlation
#' @param  w      desired confidence interval width
#' @param  p      proportion of members in one of the two subpopulations
#'
#' 
#' @references
#' \insertRef{Bonett2020a}{statpsych}
#'
#'
#' @return 
#' Returns the required sample size
#' 
#' 
#' @examples
#' size.ci.pbcor(.05, .40, .25, .73)
#'
#' # Should return:
#' # Sample size
#' #         168
#'  
#' 
#' @importFrom stats qnorm
#' @export  
size.ci.pbcor <- function(alpha, cor, w, p) {
 if (cor > .999 | cor < -.999) {stop("correlation must be between -.999 and .999")}
 z <- qnorm(1 - alpha/2)
 n <- ceiling(4*((1 - cor^2)^2)*(1 - 1.5*cor^2 + cor^2/(4*p*(1 - p)))*(z/w)^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- c("Sample size")
 rownames(out) <- ""
 return(out)
}


#  size.ci.rsqr ==============================================================
#' Sample size for a squared multiple correlation confidence interval
#'                       
#'
#' @description
#' Computes the sample size required to estimate a population squared multiple 
#' correlation in a random-x regression model with desired confidence interval
#' precision. Set the planning value of the squared multiple correlation to 1/3
#' for a conservatively large sample size. 
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
#' size.ci.rsqr(.05, .25, 5, .2)
#' 
#' # Should return:
#' # Sample size
#' #         214
#'  
#' 
#' @importFrom stats qnorm
#' @export  
size.ci.rsqr <- function(alpha, r2, s, w) {
 if (r2 > .999 | r2 < .001) {stop("squared multiple correlation must be between .001 and .999")}
 z <- qnorm(1 - alpha/2)
 n1 <- ceiling(16*(r2*(1 - r2)^2)*(z/w)^2 + s + 2)
 ci <- ci.rsqr(alpha, r2, s, n1)
 ll <- ci[1,4]                           
 ul <- ci[1,5]
 n2 <- ceiling(n1*((ul - ll)/w)^2)
 if (n2 < s + 2) {n2 = s + 2}
 ci <- ci.rsqr(alpha, r2, s, n2)
 ll <- ci[1,4]                          
 ul <- ci[1,5]
 n <- ceiling(n2*((ul - ll)/w)^2)
 if (n < s + 2) {n = s + 2}
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 rownames(out) <- ""
 return(out)
}


#  size.ci.condmean ==========================================================
#' Sample size for a conditional mean confidence interval  
#'
#'
#' @description
#' Computes the total sample size required to estimate a population conditional
#' mean of y at x = x* in a fixed-x linear regression model with desired 
#' confidence interval precision. The total sample size would be allocated to
#' the levels of the quantitative factor, and it might be necessary to increase 
#' the total sample size to give the desired sample size at each level of
#' the fixed factor. Set the error variance planning value to the largest value 
#' within a plausible range for a conservatively large sample size.
#'
#'  
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  evar   planning value of within group (error) variance
#' @param  xvar   variance of fixed predictor variable 
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
#' # Total sample size
#' #               210
#'  
#' 
#' @importFrom stats qnorm
#' @export  
size.ci.condmean <- function(alpha, evar, xvar, diff, w) {
 z <- qnorm(1 - alpha/2)
 n <- ceiling(4*(evar*(1 + diff^2/xvar))*(z/w)^2 + 1 + z^2/2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Total sample size"
 rownames(out) <- ""
 return(out)
}


#  size.ci.lc.ancova =========================================================
#' Sample size for a linear contrast confidence interval in an ANCOVA  
#'
#'
#' @description
#' Computes the sample size for each group (assuming equal sample sizes) 
#' required to estimate a population linear contrast of means in an ANCOVA model
#' with desired confidence interval precision. In a nonexperimental design,
#' the sample size is affected by the magnitude of covariate mean differences 
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
#' Returns the required sample size for each group
#' 
#' 
#' @examples
#' v <- c(1, -1)
#' size.ci.lc.ancova(.05, 1.37, 1, 0, 1.5, v)
#'
#' # Should return:
#' # Sample size per group
#' #                    21
#'  
#' 
#' @importFrom stats qnorm
#' @export  
size.ci.lc.ancova <- function(alpha, evar, s, d, w, v) {
 z <- qnorm(1 - alpha/2)
 m <- length(v) - sum(v == 0)
 n <- ceiling(4*evar*(1 + d^2/4)*(t(v)%*%v)*(z/w)^2 + s + z^2/(2*m))
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size per group"
 rownames(out) <- ""
 return(out)
}


#  size.ci.indirect ===========================================================
#' Sample size for an indirect effect confidence interval 
#'
#'
#' @description
#' Computes the approximate sample size required to estimate a population 
#' standardized indirect effect in a simple mediation model. The direct effect
#' of the independent (exogenous) variable on the response variable, controlling
#' for the mediator variable, is assumed to be negligible. 
#'
#'  
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  cor1   planning value of correlation between the independent and mediator variables
#' @param  cor2   planning value of correlation between the mediator and response variables 
#' @param  w      desired confidence interval width
#'
#' 
#' @return 
#' Returns the required sample size
#' 
#' 
#' @examples
#' size.ci.indirect(.05, .4, .5, .2)
#'
#' # Should return:
#' # Sample size
#' #         106
#'  
#' 
#' @importFrom stats qnorm
#' @export  
size.ci.indirect <- function(alpha, cor1, cor2, w) {
 if (cor1 > .999 | cor1 < -.999) {stop("cor1 must be between -.999 and .999")}
 if (cor2 > .999 | cor2 < -.999) {stop("cor2 must be between -.999 and .999")}
 z <- qnorm(1 - alpha/2)
 n <- ceiling(4*(cor2^2*(1 - cor1^2)^2 + cor1^2*(1 - cor2^2)^2)*(z/w)^2 + 3)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 rownames(out) <- ""
 return(out)
}


#  size.ci.cronbach2 =======================================================
#' Sample size for a 2-group Cronbach reliability difference confidence
#' interval
#'
#'
#' Computes the sample size per group (assuming equal sample sizes) required
#' to estimate a difference in population Cronbach reliability coefficients
#' with desired precision in a 2-group design. 
#'
#'
#' @param  alpha  alpha level for hypothesis test 
#' @param  rel1   reliability planning value for group 1
#' @param  rel2   reliability planning value for group 2
#' @param  r      number of measurements (items, raters, forms)
#' @param  w      desired confidence interval width
#'
#'
#' @return 
#' Returns the required sample size for each group
#'
#'
#' @references
#' \insertRef{Bonett2015}{statpsych}
#'
#'
#' @examples
#' size.ci.cronbach2(.05, .85, .70, 8, .15)
#'
#' # Should return:
#' # Sample size per group
#' #                   180
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.ci.cronbach2 <- function(alpha, rel1, rel2, r, w) {
 if (rel1 > .999 | rel1 < .001) {stop("rel1 must be between .001 and .999")}
 if (rel2 > .999 | rel2 < .001) {stop("rel2 must be between .001 and .999")}
 z <- qnorm(1 - alpha/2)
 n0 <- ceiling((8*r/(r - 1))*((1 - rel1)^2 + (1 - rel2)^2)*(z/w)^2 + 2)
 b <- log(n0/(n0 - 1))
 LL1 <- 1 - exp(log(1 - rel1) - b + z*sqrt(2*r/((r - 1)*(n0 - 2))))
 UL1 <- 1 - exp(log(1 - rel1) - b - z*sqrt(2*r/((r - 1)*(n0 - 2))))
 LL2 <- 1 - exp(log(1 - rel2) - b + z*sqrt(2*r/((r - 1)*(n0 - 2))))
 UL2 <- 1 - exp(log(1 - rel2) - b - z*sqrt(2*r/((r - 1)*(n0 - 2))))
 ll <- rel1 - rel2 - sqrt((rel1 - LL1)^2 + (UL2 - rel2)^2)
 ul <- rel1 - rel2 + sqrt((UL1 - rel1)^2 + (rel2 - LL2)^2)
 w0 <- ul - ll
 n <- ceiling((n0 - 2)*(w0/w)^2 + 2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size per group"
 rownames(out) <- ""
 return(out)
}


#  size.ci.mape ==============================================================
#' Sample size for a mean absolute prediction error confidence interval
#'
#'
#' Computes the sample size required to estimate a population mean absolute 
#' prediction error for a general linear model with desired confidence interval
#' precision. Setting s = 0 gives the sample size requirement for a mean absolute 
#' deviation in a one-group design. This function assumes that the prediction
#' errors have an approximate normal distribution.
#'
#'
#' @param  alpha  alpha value for 1-alpha confidence 
#' @param  mape   mean absolute prediction error planning value
#' @param  s      number of predictor variables
#' @param  w      desired confidence interval width
#'
#'
#' @return 
#' Returns the required sample size
#'
#'
#' @examples
#' size.ci.mape(.05, 4.5, 5, 2)
#'
#' # Should return:
#' # Sample size
#' #          57
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.ci.mape <- function(alpha, mape, s, w) {
 z <- qnorm(1 - alpha/2)
 n0 <- ceiling(2.28*mape^2*(z/w)^2) + s
 df <- n0 - s - 1
 c <- n0/(n0 - (s + 2)/2)
 ll <- exp(log(c*mape) - z*sqrt(.57/df))
 ul <- exp(log(c*mape) + z*sqrt(.57/df))
 w0 <- ul - ll
 n1 <- ceiling((n0 - s)*(w0/w)^2 + s)
 df <- n1 - s - 1
 c <- n1/(n1 - (s + 2)/2)
 ll <- exp(log(c*mape) - z*sqrt(.57/df))
 ul <- exp(log(c*mape) + z*sqrt(.57/df))
 w0 <- ul - ll
 n <- ceiling((n1 - s)*(w0/w)^2 + s)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 rownames(out) <- ""
 return(out)
}


#  size.ci.cor2 ==============================================================
#' Sample size for a 2-group Pearson correlation difference confidence 
#' interval 
#'
#'                     
#' @description
#' Computes the sample size required to estimate a difference in population 
#' Pearson or partial correlations with desired confidence interval precision
#' in a 2-group design. Set the correlation planning values to the smallest
#' absolute values within their plausible ranges for a conservatively large
#' sample size.
#'
#'  
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  cor1   correlation planning value for group 1
#' @param  cor2   correlation planning value for group 2
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
#' size.ci.cor2(.05, .8, .5, .2)
#'
#' # Should return:
#' # Sample size per group
#' #                   271
#'  
#' 
#' @importFrom stats qnorm
#' @export  
size.ci.cor2 <- function(alpha, cor1, cor2, w) {
 if (cor1 > .999 | cor1 < -.999) {stop("correlation must be between -.999 and .999")}
 if (cor2 > .999 | cor2 < -.999) {stop("correlation must be between -.999 and .999")}
 z <- qnorm(1 - alpha/2)
 n1 <- ceiling(4*((1 - cor1^2)^2 + (1 - cor2^2)^2)*(z/w)^2 + 3)
 ci <- ci.cor2(alpha, cor1, cor2, n1, n1)
 ll <- ci[1,3]                                  
 ul <- ci[1,4]
 n2 <- ceiling(n1*((ul - ll)/w)^2)
 ci <- ci.cor2(alpha, cor1, cor2, n2, n2)
 ll <- ci[1,3]                                  
 ul <- ci[1,4]
 n <- ceiling(n2*((ul - ll)/w)^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size per group"
 rownames(out) <- ""
 return(out)
}


# size.ci.spear2 =============================================================
#' Sample size for a 2-group Spearman correlation difference confidence 
#' interval 
#'
#'                     
#' @description
#' Computes the sample size required to estimate a difference in population 
#' Spearman correlations with desired confidence interval precision in a 2-group 
#' design. Set the correlation planning values to the smallest absolute values
#' within their plausible ranges for a conservatively large sample size.
#'
#'  
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  cor1   Spearman correlation planning value for group 1
#' @param  cor2   Spearman correlation planning value for group 2
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
#' size.ci.spear2(.05, .8, .5, .2)
#'
#' # Should return:
#' # Sample size per group
#' #                   314
#'  
#' 
#' @importFrom stats qnorm
#' @export  
size.ci.spear2 <- function(alpha, cor1, cor2, w) {
 if (cor1 > .999 | cor1 < -.999) {stop("correlation must be between -.999 and .999")}
 if (cor2 > .999 | cor2 < -.999) {stop("correlation must be between -.999 and .999")}
 z <- qnorm(1 - alpha/2)
 n1 <- ceiling(4*((1 + cor1^2/2)*(1 - cor1^2)^2 + (1 + cor2^2/2)*(1 - cor2^2)^2)*(z/w)^2 + 3)
 ci <- ci.spear2(alpha, cor1, cor2, n1, n1)
 ll <- ci[1,3]                                  
 ul <- ci[1,4]
 n2 <- ceiling(n1*((ul - ll)/w)^2)
 ci <- ci.spear2(alpha, cor1, cor2, n2, n2)
 ll <- ci[1,3]                                  
 ul <- ci[1,4]
 n <- ceiling(n2*((ul - ll)/w)^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size per group"
 rownames(out) <- ""
 return(out)
}


#  size.ci.cor.prior ==========================================================
#' Sample size for a Pearson correlation confidence interval using a 
#' planning value from a prior study
#'
#'                
#' @description
#' Computes the sample size required to estimate a Pearson correlation with
#' desired confidence interval precision in applications where an estimated
#' Pearson correlation from a prior study is available. The actual confidence
#' interval width in the planned study will depend on the value of the
#' estimated correlation in the planned study. An estimated correlation from
#' a prior study can be used to compute a prediction interval for the value of
#' the estimated correlation in the planned study. If the prediction interval
#' includes 0, then the correlation planning value is set to 0; otherwise, the
#' correlation planning value is set to the lower prediction limit (if the prior
#' correlation is positive) or the upper prediction limit (if the prior correlation
#' is negative). Using a larger confidence level (1 - alpha2) for the prediction 
#' interval will increase the probability that the width of the confidence interval
#' in the planned study will be less than or equal to the desired width.
#'
#' This sample size approach assumes that the population Pearson correlation 
#' that was estimated in the prior study is very similar to the population Pearson
#' correlation that will be estimated in the planned study. However, this type of
#' prior information is typically not available and the researcher must use expert
#' opinion to guess the value of the Pearson correlation that will be observed in the 
#' planned study. The \link[statpsych]{size.ci.cor} function uses a 
#' correlation planning value that is based on expert opinion regarding the 
#' likely value of the correlation estimate that will be observed in the 
#' planned study.
#'
#'
#' @param  alpha1  alpha level for 1-alpha1 confidence in the planned study
#' @param  alpha2  alpha level for the 1-alpha2 prediction interval 
#' @param  cor0    estimated correlation in prior study
#' @param  n0      sample size in prior study
#' @param  w       desired confidence interval width
#'
#'
#' @return
#' Returns the required sample size
#'
#'
#' @examples
#' size.ci.cor.prior(.05, .10, .438, 100, .2)
#'
#' # Should return:
#' # Sample size
#' #         331
#'
#'
#' @importFrom stats qnorm
#' @export                 
size.ci.cor.prior <- function(alpha1, alpha2, cor0, n0, w) {
 if (cor0 > .999 | cor0 < -.999) {stop("correlation must be between -.999 and .999")}
 ci <- ci.cor(alpha2, cor0, 0, n0)
 ll0 <- ci[1,3]                                  
 ul0 <- ci[1,4]  
 if (ll0 < 0 & ul0 > 0) {
   cor = 0
   n <- size.ci.cor(alpha1, cor, 0, w)
 } else {
   if (abs(ll0) < abs(ul0)) {cor = ll0}
   if (abs(ll0) > abs(ul0)) {cor = ul0}
   n <- size.ci.cor(alpha1, cor, 0, w)
   pi <- pi.cor(alpha2, cor0, n0, n, type = 1)
   ll <- pi[1,1]                                  
   ul <- pi[1,2]
   if (ll < 0 & ul > 0) {
     cor = 0
     n <- size.ci.cor(alpha1, cor, 0, w)
   } else {
     if (abs(ll) < abs(ul)) {cor = ll}
     if (abs(ll) > abs(ul)) {cor = ul}
     n <- size.ci.cor(alpha1, cor, 0, w)
	 pi <- pi.cor(alpha2, cor0, n0, n, type = 1)
     ll <- pi[1,1]                                  
     ul <- pi[1,2]
	 if (abs(ll) < abs(ul)) {cor = ll}
     if (abs(ll) > abs(ul)) {cor = ul}
     n <- size.ci.cor(alpha1, cor, 0, w)
   }
 }	 
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 rownames(out) <- ""
 return(out)
}


#  size.ci.ancova2 =========================================================
#' Sample size for a 2-group ANCOVA confidence interval
#'
#'
#' @description
#' Computes the sample size for each group required to estimate a mean 
#' difference in a 2-group ANCOVA model with desired confidence interval 
#' precision. In a nonexperimental design, the sample size is affected by 
#' the magnitude of covariate mean differences across groups. The covariate
#' mean differences can be approximated by specifying the largest 
#' standardized covariate mean difference of all covariates. In an 
#' experiment, this standardized mean difference should be set to 0. Set 
#' the error variance planning value to the largest value within a 
#' plausible range for a conservatively large sample size.
#'
#'  
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  evar   planning value of within group (error) variance
#' @param  s      number of covariates 
#' @param  d      largest standardized mean difference of all covariates
#' @param  w      desired confidence interval width
#' @param  R      ratio of n2/n1
#'
#' 
#' @return 
#' Returns the required sample size for each group
#' 
#' 
#' @examples
#' size.ci.ancova2(.05, 1.37, 1, 0, 1.5, 1)
#'
#' # Should return:
#' #  n1 n2
#' #  21 21
#'
#' size.ci.ancova2(.05, 1.37, 1, 0, 1.5, 2)
#'
#' # Should return:
#' #  n1 n2
#' #  16 32
#'
#' size.ci.ancova2(.05, 1.37, 1, .75, 1.5, 1)
#'
#' # Should return:
#' #  n1 n2
#' #  24 24
#'  
#' 
#' @importFrom stats qnorm
#' @export  
size.ci.ancova2 <- function(alpha, evar, s, d, w, R) {
 z <- qnorm(1 - alpha/2)
 n1 <- ceiling(4*evar*(1 + d^2/4)*(1 + 1/R)*(z/w)^2 + s + z^2/4)
 n2 <- ceiling(R*n1)
 out <- t(c(n1, n2))
 colnames(out) <- c("n1", "n2")
 rownames(out) <- ""
 return(out)
}


#  size.ci.slope.gen ==========================================================
#' Sample size for a slope confidence interval in a general statistical model  
#'
#'
#' @description
#' Computes the sample size required to estimate a slope coefficient with
#' desired confidence interval precision in any type of statistical model. 
#' This function requires a standard error estimate for the slope of interest 
#' from a prior or pilot study and the sample size that was used in the prior 
#' or pilot study. This function can be used for both unstandardized and
#' standardized slopes. This function also can be used for both unstandardized 
#' and standardized factor loadings in a confirmatory factor analysis model.
#' This function will soon be replaced with size.ci.gen.
#'
#'  
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  se     standard error of slope from prior/pilot study
#' @param  n0     sample size used in prior/pilot study 
#' @param  w      desired confidence interval width
#'
#' 
#' @return 
#' Returns the required sample size
#' 
#' 
#' @examples
#' size.ci.slope.gen(.05, 3.15, 50, 5)
#'
#' # Should return:
#' #  Sample size
#' #          305
#'  
#' 
#' @importFrom stats qnorm
#' @export  
size.ci.slope.gen <- function(alpha, se, n0, w) {
 z <- qnorm(1 - alpha/2)
 n <- ceiling(4*n0*se^2*(z/w)^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 rownames(out) <- ""
 return(out)
}


#  size.ci.gen ============================================================
#' Sample size for a confidence interval for any type of parameter
#'
#'
#' @description
#' Computes the sample size required to estimate a single population parameter
#' with desired precision using a standard error for the parameter estimate 
#' from a prior or pilot study. This function can be used with any type of 
#' parameter where the standard error of the parameter estimate is a function
#' of the square root of the sample size (most parameter estimates have this
#' property). This function also assumes that the sampling distribution of the 
#' parameter estimate is approximately normal in large samples.
#'
#'
#' @param  alpha  alpha level for hypothesis test 
#' @param  se     standard error of parameter estimate from prior/pilot study
#' @param  n0     sample size of prior/pilot study
#' @param  w      desired confidence interval width
#'
#'
#' @return 
#' Returns the required sample size 
#'
#'
#' @examples
#' size.ci.gen(.05, 2.89, 30, 8)
#'
#' # Should return:
#' # Sample size
#' #          61
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.ci.gen <- function(alpha, se, n0, w) {
 if (w <= 0) {stop("width must be a positive value")}
 warn <- "Warning: alpha level is typically less than .25"
 if (alpha > .25) {message(warn)}
 za <- qnorm(1 - alpha/2)
 n <- ceiling(4*n0*se^2*(za/w)^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 rownames(out) <- ""
 return(out)
}


#  size.ci.gen2 ============================================================
#' Sample size for a confidence interval for the difference of any type of 
#' parameter
#'
#'
#' @description
#' Computes the sample size required to estimate a difference in population
#' parameters with desired precision in a 2-group design using a standard
#' error for a parameter estimate from a prior or pilot study. This function
#' can be used with any type of parameter where the standard error of the 
#' parameter estimate is a function of the square root of the sample size 
#' (most parameter estimates have this property). This function also assumes 
#' that the sampling distribution of the parameter estimate is approximately
#' normal in large samples. Set R = 1 for equal sample sizes.
#'
#'
#' @param  alpha  alpha level for hypothesis test 
#' @param  se     standard error of parameter estimate from prior/pilot study
#' @param  n0     sample size of prior/pilot study
#' @param  w      desired confidence interval width
#' @param  R      n2/n1 ratio
#'
#'
#' @return 
#' Returns the required sample size for each group 
#'
#'
#' @examples
#' size.ci.gen2(.05, .175, 30, .8, 1)
#'
#' # Should return:
#' # n1  n2
#' # 45  45
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.ci.gen2 <- function(alpha, se, n0, w, R) {
 if (w <= 0) {stop("width must be a positive value")}
 warn <- "Warning: alpha level is typically less than .25"
 if (alpha > .25) {message(warn)}
 za <- qnorm(1 - alpha/2)
 n1 <- ceiling(4*(1 + 1/R)*n0*se^2*(za/w)^2)
 n2 <- ceiling(R*n1)
 out <- t(c(n1, n2))
 colnames(out) <- c("n1", "n2")
 rownames(out) <- ""
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
#' @param  h       null hypothesis value of slope  
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
#' # Total sample size
#' #               100
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
 colnames(out) <- "Total sample size"
 rownames(out) <- ""
 return(out)
}


#  size.test.cor ==============================================================
#' Sample size for a test of a Pearson or partial correlation 
#'
#'
#' @description
#' Computes the sample size required to test a population Pearson or a partial
#' correlation with desired power. Set s = 0 for a Pearson correlation. 
#'
#'  
#' @param  alpha   alpha level for hypothesis test
#' @param  pow     desired power
#' @param  cor     planning value of correlation
#' @param  s       number of control variables
#' @param  h       null hypothesis value of correlation  
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
#' # Sample size
#' #          48
#'  
#' 
#' @importFrom stats qnorm
#' @export  
size.test.cor <- function(alpha, pow, cor, s, h) {
 if (cor > .999 | cor < -.999) {stop("correlation must be between -.999 and .999")}
 za <- qnorm(1 - alpha/2)
 zb <- qnorm(pow)
 zr <- log((1 + cor)/(1 - cor))/2
 zo <- log((1 + h)/(1 - h))/2
 es <- zr - zo
 n <- ceiling((za + zb)^2/es^2 + s + 3)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 rownames(out) <- ""
 return(out)
}


#  size.interval.cor =========================================================
#' Sample size for an interval test of a Pearson or partial correlation 
#'
#'
#' @description
#' Computes the sample size required to perform an interval test for a 
#' population Pearson or a partial correlation with desired power where the 
#' interval midpoint is equal to zero. This function can be used to plan a study
#' where the goal is to show that the population correlation is small. Set s = 0
#' for a Pearson correlation. The correlation planning value must be a value 
#' within the hypothesized interval. 
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
#' # Sample size
#' #         360
#'  
#' 
#' @importFrom stats qnorm
#' @export  
size.interval.cor <- function(alpha, pow, cor, s, h) {
 if (cor > .999 | cor < -.999) {stop("correlation must be between -.999 and .999")}
 if (h <= abs(cor)) {stop("correlation must be between -h and h")}
 za <- qnorm(1 - alpha)
 zb <- qnorm(1 - (1 - pow)/2)
 zr <- log((1 + cor)/(1 - cor))/2
 zh <- log((1 + h)/(1 - h))/2
 es <- zh - zr
 n <- ceiling((za + zb)^2/es^2 + s + 3)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 rownames(out) <- ""
 return(out)
}


#  size.test.cor2 ==============================================================
#' Sample size for a test of equal Pearson or partial correlation in a 2-group
#' design
#'
#'
#' @description
#' Computes the sample size required to test equality of two Pearson or partial
#' correlation with desired power in a 2-group design. Set s = 0 for a Pearson 
#' correlation. Set R = 1 for equal sample sizes.
#'
#'  
#' @param  alpha   alpha level for hypothesis test
#' @param  pow     desired power
#' @param  cor1    correlation planning value for group 1
#' @param  cor2    correlation planning value for group 2
#' @param  s       number of control variables
#' @param  R       n2/n1 ratio
#'
#' 
#' @return 
#' Returns the required sample size for each group
#' 
#' 
#' @examples
#' size.test.cor2(.05, .8, .4, .2, 0, 1)
#'
#' # Should return:
#' #  n1  n2
#' # 325 325
#'
#' size.test.cor2(.05, .8, .4, .2, 0, 2)
#'
#' # Should return:
#' #  n1  n2
#' # 245 490
#'  
#' 
#' @importFrom stats qnorm
#' @export  
size.test.cor2 <- function(alpha, pow, cor1, cor2, s, R) {
 if (cor1 > .999 | cor1 < -.999) {stop("cor1 must be between -.999 and .999")}
 if (cor2 > .999 | cor2 < -.999) {stop("cor2 must be between -.999 and .999")}
 za <- qnorm(1 - alpha/2)
 zb <- qnorm(pow)
 zr1 <- log((1 + cor1)/(1 - cor1))/2
 zr2 <- log((1 + cor2)/(1 - cor2))/2
 es <- zr1 - zr2
 n1 <- ceiling((1 + 1/R)*(za + zb)^2/es^2 + s + 3)
 n2 <- ceiling(R*n1)
 out <- t(c(n1, n2))
 colnames(out) <- c("n1", "n2")
 rownames(out) <- ""
 return(out)
}


#  size.test.lc.ancova ==========================================================
#' Sample size for a mean linear contrast test in an ANCOVA 
#'
#'
#' @description
#' Computes the sample size for each group (assuming equal sample sizes) required
#' to test a linear contrast of population means in an ANCOVA model with desired
#' power. In a nonexperimental design, the sample size is affected by the magnitude
#' of covariate mean differences across groups. The covariate mean differences can  
#' be approximated by specifying the largest standardized covariate mean difference 
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
#' Returns the required sample size for each group
#' 
#' 
#' @examples
#' v <- c(.5, .5, -1)
#' size.test.lc.ancova(.05, .9, 1.37, .7, 1, 0, v)
#'
#' # Should return:
#' # Sample size per group
#' #                    46
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.test.lc.ancova <- function(alpha, pow, evar, es, s, d, v) {
 za <- qnorm(1 - alpha/2)
 zb <- qnorm(pow)
 m <- length(v) - sum(v == 0)
 n <- ceiling((evar*(1 + d^2/4)*t(v)%*%v)*(za + zb)^2/es^2 + s + za^2/(2*m))
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size per group"
 rownames(out) <- ""
 return(out)
}


#  size.test.cronbach2 ========================================================
#' Sample size to test equality of Cronbach reliability coefficients in a 
#' 2-group design
#'
#'
#' @description
#' Computes the sample size required to test a difference in population
#' Cronbach reliability coefficients with desired power in a 2-group design. 
#'
#'
#' @param  alpha  alpha level for hypothesis test 
#' @param  pow    desired power
#' @param  rel1   reliability planning value for group 1
#' @param  rel2   reliability planning value for group 2
#' @param  r      number of measurements (items, raters, forms)
#'
#'
#' @return 
#' Returns the required sample size for each group
#'
#'
#' @references
#' \insertRef{Bonett2015}{statpsych}
#'
#'
#' @examples
#' size.test.cronbach2(.05, .80, .85, .70, 8)
#'
#' # Should return:
#' # Sample size per group
#' #                    77
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.test.cronbach2 <- function(alpha, pow, rel1, rel2, r) {
 if (rel1 > .999 | rel1 < .001) {stop("rel1 must be between .001 and .999")}
 if (rel2 > .999 | rel2 < .001) {stop("rel2 must be between .001 and .999")}
 za <- qnorm(1 - alpha/2)
 zb <- qnorm(pow)
 e <- (1 - rel1)/(1 - rel2)
 n <- ceiling((4*r/(r - 1))*(za + zb)^2/log(e)^2 + 2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size per group"
 rownames(out) <- ""
 return(out)
}


#  size.test.ancova2 ==========================================================
#' Sample size for a 2-group ANCOVA hypothesis test
#'
#'
#' @description
#' Computes the sample size for each group required to test a mean difference
#' in an ANCOVA model with desired power in a 2-group design. In a 
#' nonexperimental design, the sample size is affected by the magnitude of 
#' covariate mean differences across groups. The covariate mean differences can
#' be approximated by specifying the largest standardized covariate mean 
#' difference across of all covariates. In an experiment, this standardized 
#' mean difference is set to 0. Set the error variance planning value to the 
#' largest value within a plausible range for a conservatively large sample 
#' size.
#'
#'  
#' @param  alpha   alpha level for hypothesis test
#' @param  pow     desired power
#' @param  evar    planning value of within-group (error) variance
#' @param  es      planning value of mean difference
#' @param  s       number of covariates 
#' @param  d       largest standardized mean difference of all covariates
#' @param  R       n2/n1 rartio
#'
#' 
#' @return 
#' Returns the required sample size for each group
#' 
#' 
#' @examples
#' size.test.ancova2(.05, .9, 1.37, .7, 1, 0, 1)
#'
#' # Should return:
#' #  n1 n2
#' #  61 61
#'
#' size.test.ancova2(.05, .9, 1.37, .7, 1, 0, 2)
#'
#' # Should return:
#' #  n1 n2
#' #  47 94
#'
#' size.test.ancova2(.05, .9, 1.37, .7, 1, .5, 1)
#'
#' # Should return:
#' #  n1 n2
#' #  65 65
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.test.ancova2 <- function(alpha, pow, evar, es, s, d, R) {
 if (es == 0) {stop("effect size cannot be zero")}
 za <- qnorm(1 - alpha/2)
 zb <- qnorm(pow)
 n1 <- ceiling((evar*(1 + d^2/4)*(1 + 1/R))*(za + zb)^2/es^2 + s + za^2/4)
 n2 <- ceiling(R*n1)
 out <- t(c(n1, n2))
 colnames(out) <- c("n1", "n2")
 rownames(out) <- ""
 return(out)
}


#  size.test.slope.gen ========================================================
#' Sample size for a slope hypothesis test in a general statistical model  
#'
#'
#' @description
#' Computes the sample size required to test a null hypothesis with desired
#' power that a population slope coefficient in any general statistical model 
#' is equal to zero. This function requires a standard error estimate for the 
#' slope of interest from a prior or pilot study and the sample size that was
#' used in the prior or pilot study. This function can be used for both 
#' unstandardized and standardized slopes. This function also can be used for 
#' both unstandardized and standardized factor loadings in a confirmatory
#' factor analysis model. This function will soon be replaced with size.test.gen.
#'
#'  
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  pow    desired power
#' @param  se     standard error of slope from prior/pilot study
#' @param  n0     sample size used in prior/pilot study 
#' @param  b      planning value of population slope
#'
#' 
#' @return 
#' Returns the required sample size
#' 
#' 
#' @examples
#' size.test.slope.gen(.05, .8, 3.15, 50, 5)
#'
#' # Should return:
#' #  Sample size
#' #          156
#'  
#' 
#' @importFrom stats qnorm
#' @export  
size.test.slope.gen <- function(alpha, pow, se, n0, b) {
 if (b == 0) {stop("slope planning value cannot be zero")}
 za <- qnorm(1 - alpha/2)
 zb <- qnorm(pow)
 n <- ceiling(n0*se^2*(za + zb)^2/b^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 rownames(out) <- ""
 return(out)
}


#  size.test.gen ============================================================
#' Sample size for a test of any type of parameter
#'
#'
#' @description
#' Computes the sample size required to test a single population parameter with 
#' desired power using a standard error for the parameter estimate from a prior
#' or pilot study. This function can be used with any type of parameter where the
#' standard error of the parameter estimate is a function of the square root of the
#' sample size (most parameter estimates have this property). This function also 
#' assumes that the sampling distribution of the parameter estimate is 
#' approximately normal in large samples.
#'
#'
#' @param  alpha  alpha level for hypothesis test 
#' @param  pow    desired power
#' @param  se     standard error of parameter estimate from prior/pilot study
#' @param  n0     sample size of prior/pilot study
#' @param  es     planning value of parameter minus null hypothesis value
#'
#'
#' @return 
#' Returns the required sample size 
#'
#'
#' @examples
#' size.test.gen(.05, .8, 2.89, 30, 5)
#'
#' # Should return:
#' # Sample size
#' #          79
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.test.gen <- function(alpha, pow, se, n0, es) {
 if (es == 0) {stop("effect size cannot equal 0")}
 warn <- "Warning: alpha level is typically less than .25"
 if (alpha > .25) {message(warn)}
 za <- qnorm(1 - alpha/2)
 zb <- qnorm(pow)
 n <- ceiling(n0*se^2*(za + zb)^2/es^2)
 out <- matrix(n, nrow = 1, ncol = 1)
 colnames(out) <- "Sample size"
 rownames(out) <- ""
 return(out)
}


#  size.test.gen2 ============================================================
#' Sample size for a test of 2-group difference for any type of parameter
#'
#'
#' @description
#' Computes the sample size per group required to test a difference in two
#' populatation parameters with desired power using a standard error for a 
#' single parameter estimate from a prior or pilot study. This function can be
#' used with any type of parameter where the standard error of the parameter 
#' estimate is a function of the square root of the sample size (most parameter 
#' estimates have this property). This function also assumes that the sampling
#' distribution of the parameter estimate is approximately normal in large 
#' samples. Set R = 1 for equal sample sizes.
#'
#'
#' @param  alpha  alpha level for hypothesis test 
#' @param  pow    desired power
#' @param  se     standard error of parameter estimate from prior/pilot study
#' @param  n0     sample size of prior/pilot study
#' @param  es     planning value of parameter difference
#' @param  R      n2/n1 ratio
#'
#'
#' @return 
#' Returns the required sample size for each group 
#'
#'
#' @examples
#' size.test.gen2(.05, .85, .175, 30, .5, 1)
#'
#' # Should return:
#' # n1  n2
#' # 66  66
#'  
#' 
#' @importFrom stats qnorm
#' @export
size.test.gen2 <- function(alpha, pow, se, n0, es, R) {
 if (es == 0) {stop("effect size cannot equal 0")}
 warn <- "Warning: alpha level is typically less than .25"
 if (alpha > .25) {message(warn)}
 za <- qnorm(1 - alpha/2)
 zb <- qnorm(pow)
 n1 <- ceiling((1 + 1/R)*n0*se^2*(za + zb)^2/es^2)
 n2 <- ceiling(R*n1)
 out <- t(c(n1, n2))
 colnames(out) <- c("n1", "n2")
 rownames(out) <- ""
 return(out)
}


#  ======================= Power for Planned Sample Size ======================
#  power.cor =================================================================
#' Approximates the power of a correlation test for a planned sample size
#'
#'
#' @description
#' Computes the approximate power of a test for a population Pearson or partial
#' correlation test for a planned sample size. Set s = 0 for a Pearson correlation. 
#'
#'
#' @param  alpha  alpha level for hypothesis test 
#' @param  n      planned sample size
#' @param  cor    planning value of correlation 
#' @param  h      null hypothesis value of correlation 
#' @param  s      number of control variables
#'
#'
#' @return
#' Returns the approximate power of the test
#'
#'
#' @examples
#' power.cor(.05, 80, .3, 0, 0)
#'
#' # Should return:
#' #     Power
#' # 0.7751947
#'
#'
#' @importFrom stats qnorm
#' @importFrom stats pnorm
#' @export
power.cor <- function(alpha, n, cor, h, s) {
 if (cor > .999 | cor < -.999) {stop("correlation must be between -.999 and .999")}
 if (h > .999 | h < -.999) {stop("h must be between -.999 and .999")}
 za <- qnorm(1 - alpha/2)
 f1 <- log((1 + cor)/(1 - cor))/2 
 f2 <- log((1 + h)/(1 - h))/2
 z1 <- abs(f1 - f2)*sqrt(n - 3 - s) - za
 z2 <- abs(f1 - f2)*sqrt(n - 3 - s) + za
 pow1 <- pnorm(z1)
 pow2 <- 1 - pnorm(z2)
 pow <- pow1 + pow2
 out <- matrix(pow, nrow = 1, ncol = 1)
 colnames(out) <- "Power"
 rownames(out) <- ""
 return(out)
}


#  power.cor2 =================================================================
#' Approximates the power of a test for equal correlations in a 2-group design
#' for planned sample sizes
#'
#'
#' @description
#' Computes the approximate power of a test for equal population Pearson or
#' partial correlations in a 2-group design for planned sample sizes. Set
#' s = 0 for Pearson correlations. 
#'
#'
#' @param  alpha  alpha level for hypothesis test 
#' @param  n1     planned sample size for group 1
#' @param  n2     planned sample size for group 2
#' @param  cor1   correlation planning value for group 1 
#' @param  cor2   correlation planning value for group 2 
#' @param  s      number of control variables
#'
#'
#' @return
#' Returns the approximate power of the test
#'
#'
#' @examples
#' power.cor2(.05, 200, 200, .4, .2, 0)
#'
#' # Should return:
#' #     Power
#' # 0.5919682
#'
#'
#' @importFrom stats qnorm
#' @importFrom stats pnorm
#' @export
power.cor2 <- function(alpha, n1, n2, cor1, cor2, s) {
 if (cor1 > .999 | cor1 < -.999) {stop("cor1 must be between -.999 and .999")}
 if (cor2 > .999 | cor2 < -.999) {stop("cor2 must be between -.999 and .999")}
 za <- qnorm(1 - alpha/2)
 f1 <- log((1 + cor1)/(1 - cor1))/2
 f2 <- log((1 + cor2)/(1 - cor2))/2
 z1 <- abs(f1 - f2)/sqrt(1/(n1 - 3 - s) + 1/(n2 - 3 - s)) - za
 z2 <- abs(f1 - f2)/sqrt(1/(n1 - 3 - s) + 1/(n2 - 3 - s)) + za
 pow1 <- pnorm(z1)
 pow2 <- 1 - pnorm(z2)
 pow <- pow1 + pow2
 out <- matrix(pow, nrow = 1, ncol = 1)
 colnames(out) <- "Power"
 rownames(out) <- ""
 return(out)
}


# ============================= Miscellaneous =================================
#  slope.contrast =============================================================
#' Contrast coefficients for the slope of a quantitative factor
#'
#'
#' @description
#' Computes the contrast coefficients that are needed to estimate the slope of
#' a line in a one-factor design with a quantitative factor.
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
 colnames(out) <- "Coefficient"
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
 if (cor > .999 | cor < -.999) {stop("correlation must be between -.999 and .999")}
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
#' * 3 = leptokurtic (skewness = 0 and excess kurtosis = 6)
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
 if (cor > .999 | cor < -.999) {stop("correlation must be between -.999 and .999")}
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
#' * 3 = leptokurtic (skewness = 0 and excess kurtosis = 6)
#' * 4 = moderate skew (skewness = 1 and excess kurtosis = 1.5)
#' * 5 = large skew (skewness = 2 and excess kurtosis = 6)
#' @param   rep      number of Monte Carlo samples
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
 if (cor > .999 | cor < -.999) {stop("correlation must be between -.999 and .999")}
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


#  adj.se =====================================================================
#' Adjusted standard errors for slope coefficients in an exploratory analysis
#'                              
#'
#' @description
#' Computes adjusted standard errors in a general linear model after one or 
#' more predictor variables with nonsignificant slopes have been dropped from 
#' the model. The adjusted standard errors are then used to compute adjusted 
#' t-values, p-values, and confidence intervals. The mean square error and
#' error degrees of freedom from the full model are used to compute the 
#' adjusted standard errors. These adjusted results are less susceptible to
#' the negative effects of an exploratory model selection. 
#'
#'  
#' @param  alpha  alpha level for 1-alpha confidence
#' @param  mse1   mean squared error in full model
#' @param  mse2   mean squared error in selected model
#' @param  dfe1   error df in full model
#' @param  se     vector of slope standard errors in selected model
#' @param  b      vector of estimated slopes in selected model
#'
#' 
#' @return 
#' Returns adjusted standard error, t-statistic, p-value, and confidence interval
#' for each slope coefficient
#' 
#' 
#' @examples
#' se <- c(1.57, 3.15, 0.982)
#' b <- c(3.78, 8.21, 2.99)
#' adj.se(.05, 10.26, 8.37, 114, se, b)
#'
#' # Should return:
#' #      Estimate   adj SE        t  df           p        LL        UL
#' # [1,]     3.78 1.738243 2.174609 114 0.031725582 0.3365531  7.223447
#' # [2,]     8.21 3.487559 2.354082 114 0.020279958 1.3011734 15.118827
#' # [3,]     2.99 1.087233 2.750102 114 0.006930554 0.8362007  5.143799
#'  
#' 
#' @importFrom stats qt
#' @importFrom stats pt
#' @export  
adj.se <- function(alpha, mse1, mse2, dfe1, se, b) {
 s <- length(b)
 df <- rep(dfe1, s)
 tcrit <- qt(1 - alpha/2, dfe1)
 adjse <- se*sqrt(mse1/mse2)
 t <- b/adjse
 p <- 2*(1 - pt(abs(t),dfe1))
 ll <- b - tcrit*adjse
 ul <- b + tcrit*adjse
 out <- matrix(c(t(b), t(adjse), t(t), t(df), t(p), t(ll), t(ul)), s, 7)
 colnames(out) <- c("Estimate", "adj SE", "t", "df", "p", "LL", "UL")
 return(out)
}


#  fitindices =================================================================
#' SEM fit indices
#'
#'
#' @description
#' Computes the normed fit index (NFI), adjusted normed fit index (adj NFI),
#' comparative fit index (CFI), Tucker-Lewis fit index (TLI), and root mean
#' square error of approximation index (RMSEA). Of the first four indices, the
#' adj NFI index is recommended because it has smaller sampling variability
#' than CFI and TLI and less negative bias than NFI.
#'
#'  
#' @param  chi1   chi-square test statistic for full model
#' @param  df1    degrees of freedom for full model
#' @param  chi2   chi-square test statistic for reduced model
#' @param  df2    degrees of freedom for reduced model
#' @param  n      sample size
#'
#' 
#' @return 
#' Returns NFI, adj NFI, CFI, TLI, and RMSEA
#' 
#' 
#' @examples
#' fitindices(14.21, 10, 258.43, 20, 300)
#'
#' # Should return:
#' #        NFI   adj NFI       CFI       TLI      RMSEA
#' #  0.9450141 0.9837093 0.9823428 0.9646857 0.03746109
#'  
#' 
#' @export  
fitindices <- function(chi1, df1, chi2, df2, n) {
 if (chi2 == 0) {stop("chi2 must be a positive value")}
 nfi <- 1 - chi1/chi2
 d1 <- chi1 - df1
 if (d1 < 0) {d1 = 0}
 adjnfi <- 1 - d1/chi2
 rmsea <- sqrt(d1/(n*df1))
 if (chi2 - df2 > 0) 
  {cfi <- 1 - d1/(chi2 - df2)}
 else
  {cfi <- 0}
 d2 <- chi2/df2 - chi1/df1
 d3 <- chi2/df2 - 1
 if (d2 < 0) {d2 = 0}
 if (d3 < 0) {d3 = 0}
 if (d3 > 0)
  {tli <- d2/d3}
 else
  {tli <- 0}
 out <- t(c(nfi, adjnfi, cfi, tli, rmsea))
 colnames(out) <- c("NFI", "adj NFI", "CFI", "TLI", "RMSEA")
 rownames(out) <- ""
 return(out)
}

