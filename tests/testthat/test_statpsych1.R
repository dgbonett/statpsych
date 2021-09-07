library(statpsych)

test_that("ci.mean1 returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "LL", "UL")
  
  res <- ci.mean1(.05, 24.5, 3.65, 40)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.stdmean1 returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "LL", "UL")
  
  res <- ci.stdmean1(.05, 24.5, 3.65, 40, 20)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.mean2 returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "t", "df", "p", "LL", "UL")
  
  res <- ci.mean2(.05, 15.4, 10.3, 2.67, 2.15, 30, 20)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(2, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.lc.mean.bs returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "t", "df", "p", "LL", "UL")
  
  m <- c(33.5, 37.9, 38.0, 44.1)
  sd <- c(3.84, 3.84, 3.65, 4.98)
  n <- c(10,10,10,10)
  c <- c(.5, .5, -.5, -.5)
  res <- ci.lc.mean.bs(.05, m, sd, n, c)

  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(2, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.tukey returns valid matrix", {
  colnames_expected <- c("", "", "diff", "SE", "t", "df", "p", "LL", "UL")

  m <- c(12.86, 17.57, 26.29)
  sd <- c(3.185, 3.995, 3.773)
  n <- c(20, 20, 20)
  res <- ci.tukey(.05, m, sd, n)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(3, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.ratio.mean2  returns valid matrix", {
  colnames_expected <- c("Mean1", "Mean2", "Mean1/Mean2", "LL", "UL")
  
  y2 <- c(32, 39, 26, 35, 43, 27, 40, 37, 34, 29, 49, 42, 40)
  y1 <- c(36, 44, 47, 42, 49, 39, 46, 31, 33, 48)
  res <- ci.ratio.mean2(.05, y1, y2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.stdmean2  returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "LL", "UL")
  
  res <- ci.stdmean2(.05, 35.1, 26.7, 7.32, 6.98, 30, 30)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(4, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.stdmean.strat returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "LL", "UL")
  
  res <- ci.stdmean.strat(.05, 30.2, 30.8, 10.5, 11.2, 200, 200, .533)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(3, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.lc.stdmean.bs returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "LL", "UL")
  
  m <- c(33.5, 37.9, 38.0, 44.1)
  sd <- c(3.84, 3.84, 3.65, 4.98)
  n <- c(10,10,10,10)
  v <- c(.5, .5, -.5, -.5)
  res <- ci.lc.stdmean.bs(.05, m, sd, n, v)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(2, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.mean.ps returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "t", "df", "p", "LL", "UL")
  
  res <- ci.mean.ps(.05, 58.2, 51.4, 7.43, 8.92, .537, 30)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.ratio.mean.ps returns valid matrix", {
  colnames_expected <- c("Mean1", "Mean2", "Mean1/Mean2", "LL", "UL")
  
  y1 <- c(3.3, 3.6, 3.0, 3.1, 3.9, 4.2, 3.5, 3.3)
  y2 <- c(3.0, 3.1, 2.7, 2.6, 3.2, 3.8, 3.2, 3.0)
  res <- ci.ratio.mean.ps(.05, y1, y2)

  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.stdmean.ps returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "LL", "UL")
  
  res <- ci.stdmean.ps(.05, 110.4, 102.1, 15.3, 14.6, .75, 25)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(3, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.lc.stdmean.ws returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "LL", "UL")
  
  m <- c(33.5, 37.9, 38.0, 44.1)
  sd <- c(3.84, 3.84, 3.65, 4.98)
  q <- c(.5, .5, -.5, -.5)
  res <- ci.lc.stdmean.ws(.05, m, sd, .672, 20, q)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})



test_that("ci.mad1 returns valid matrix", {
  colnames_expected <- c("MAD", "LL", "UL")
  
  y <- c(30, 20, 15, 10, 10, 60, 20, 25, 20, 30, 10, 5, 50, 40, 
         20, 10, 0, 20, 50)
  res <- ci.mad1(.05, y)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.ratio.mad2 returns valid matrix", {
  colnames_expected <- c("MAD1", "MAD2", "MAD1/MAD2", "LL", "UL")
  
  y1 <- c(32, 39, 26, 35, 43, 27, 40, 37, 34, 29)
  y2 <- c(36, 44, 47, 42, 49, 39, 46, 31, 33, 48)
  res <- ci.ratio.mad2(.05, y1, y2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.ratio.mad.ps returns valid matrix", {
  colnames_expected <- c("MAD1", "MAD2", "MAD1/MAD2", "LL", "UL")
  
  y2 <- c(21, 4, 9, 12, 35, 18, 10, 22, 24, 1, 6, 8, 13, 16, 19)
  y1 <- c(67, 28, 30, 28, 52, 40, 25, 37, 44, 10, 14, 20, 28, 40, 51)
  res <- ci.ratio.mad.ps(.05, y1, y2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.cod1 returns valid matrix", {
  colnames_expected <- c("COD", "LL", "UL")
  
  y <- c(30, 20, 15, 10, 10, 60, 20, 25, 20, 30, 10, 5, 50, 40,
         20, 10, 0, 20, 50)
  res <- ci.cod1(.05, y)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.median1 returns valid matrix", {
  colnames_expected <- c("Median", "SE", "LL", "UL")
  
  y <- c(30, 20, 15, 10, 10, 60, 20, 25, 20, 30, 10, 5, 50, 40,
         20, 10, 0, 20, 50)
  res <- ci.median1(.05, y)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.median2 returns valid matrix", {
  colnames_expected <- c("Median1", "Median2", "Median1-Median2", "SE", "LL", "UL")
  
  y1 <- c(32, 39, 26, 35, 43, 27, 40, 37, 34, 29)
  y2 <- c(36, 44, 47, 42, 49, 39, 46, 31, 33, 48)
  res <- ci.median2(.05, y1, y2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.ratio.median2 returns valid matrix", {
  colnames_expected <- c("Median1", "Median2", "Median1/Median2", "LL", "UL")
  
  y2 <- c(32, 39, 26, 35, 43, 27, 40, 37, 34, 29, 49, 42, 40)
  y1 <- c(36, 44, 47, 42, 49, 39, 46, 31, 33, 48)
  res <- ci.ratio.median2(.05, y1, y2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.lc.median.bs returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "LL", "UL")
  
  m <- c(46.13, 29.19, 30.32, 49.15)
  se <- c(6.361, 5.892, 4.887, 6.103)
  v <- c(1, -1, -1, 1)
  res <- ci.lc.median.bs(.05, m, se, v)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.median.ps returns valid matrix", {
  colnames_expected <- c(
    "Median1", "Median2", "Median1-Median2", "SE", "LL", "UL", "SE1", "SE2", "cov"
  )
  
  y1 <- c(21, 4, 9, 12, 35, 18, 10, 22, 24, 1, 6, 8, 13, 16, 19)
  y2 <- c(67, 28, 30, 28, 52, 40, 25, 37, 44, 10, 14, 20, 28, 40, 51)
  res <- ci.median.ps(.05, y1, y2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.ratio.median.ps returns valid matrix", {
  colnames_expected <- c(
    "Median1", "Median2", "Median1/Median2", "LL", "UL"
  )
  
  y1 <- c(21, 4, 9, 12, 35, 18, 10, 22, 24, 1, 6, 8, 13, 16, 19)
  y2 <- c(67, 28, 30, 28, 52, 40, 25, 37, 44, 10, 14, 20, 28, 40, 51)
  res <- ci.ratio.median.ps(.05, y1, y2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.mann returns valid matrix", {
  colnames_expected <- c(
   "Estimate", "SE", "LL", "UL"
  )
  
  y1 <- c(32, 39, 26, 35, 43, 27, 40, 37, 34, 29)
  y2 <- c(36, 44, 47, 42, 49, 39, 46, 31, 33, 48)
  res <- ci.mann(.05, y1, y2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.random.anova1 returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "LL", "UL"
  )
  
  m <- c(56.1, 51.2, 60.3, 68.2, 48.9, 70.5)
  sd <- c(9.45, 8.79, 9.71, 8.90, 8.31, 9.75)
  res <- ci.random.anova1(.05, m, sd, 20)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(4, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.cronbach returns valid matrix", {
  colnames_expected <- c(
    "LL", "UL"
  )
  
  res <- ci.cronbach(.05, .85, 7, 89)
  
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("size.ci.mean1 returns valid number", {
  colnames_expected <- c(
    "Sample size"
  )
  
  
  res <- size.ci.mean1(.05, 264.4, 10)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  testthat::expect_equal(res[[1,1]], 43)
})


test_that("size.ci.mean2 returns valid matrix", {
  colnames_expected <- c(
    "n1", "n2"
  )
  
  res <- size.ci.mean2(.05, 37.1, 5, 1)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("size.ci.stdmean2 returns valid matrix", {
  colnames_expected <- c(
    "n1", "n2"
  )
  
  res <- size.ci.stdmean2(.05, .75, .5, 1)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("size.ci.ratio.mean2 returns valid matrix", {
  colnames_expected <- c(
    "n1", "n2"
  )
  
  res <- size.ci.ratio.mean2(.05, .4, 3.5, 3.1, 1.2, 2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("size.ci.lc.mean.bs returns valid number", {

  v <- c(.5, .5, -1)
  res <- size.ci.lc.mean.bs(.05, 5.62, 2.0, v)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1]], 34)
})


test_that("size.ci.stdmean.ps returns valid number", {
  res <- size.ci.mean.ps(.05, 265, .8, 10)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1]], 19)
})


test_that("size.ci.ratio.mean2 returns valid number", {
  res <- size.ci.stdmean.ps(.05, 1, .65, .6)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1]], 46)
})


test_that("size.ci.ratio.mean.ps returns valid number", {
  res <- size.ci.ratio.mean.ps(.05, 400, 150, 100, .7, 1.2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1]], 21)
})


test_that("size.ci.lc.stdmean.ws returns valid matrix", {
  q <- c(.5, .5, -.5, -.5)
  res <- size.ci.lc.mean.ws(.05, 265, .8, 10, q)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 11)
})


test_that("size.ci.lc.mean.ws returns valid matrix", {
  q <- c(.5, .5, -.5, -.5)
  res <- size.ci.lc.stdmean.ws(.05, 1, .7, .6, q)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 26)
})


test_that("size.ci.cronbach returns valid number", {
  res <- size.ci.cronbach(.05, .85, 5, .1)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 89)
})


test_that("size.ci.second returns valid number", {
  res <- size.ci.second(20, 5.3, 2.5)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 70)
})


test_that("size.test.mean1 returns valid number", {
  res <- size.test.mean1(.05, .9, 80.5, 7)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 20)
})


test_that("size.test.mean2 returns valid matrix", {
  colnames_expected <- c(
    "n1", "n2"
  )
  
  res <- size.test.mean2(.05, .95, 100, 10, 1) 
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("size.test.lc.mean.bs returns valid matrix", {
  v <- c(1, -1, -1, 1)
  res <- size.test.lc.mean.bs(.05, .90, 27.5, 5, v)

  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 47)  
    
})


test_that("size.equiv.mean2 returns valid matrix", {

  res <- size.equiv.mean2(.10, .80, 15, 2, 4)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 50)
})


test_that("size.supinf.mean2 returns valid matrix", {
  res <- size.supinf.mean2(.05, .80, 225, 9, 4)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 143)
})


test_that("size.test.mean.ps returns valid number", {
  res <- size.test.mean.ps(.05, .80, 1.25, .5, .75) 
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 22)
})


test_that("size.test.lc.mean.ws returns valid matrix", {
  q <- c(.5, .5, -.5, -.5)
  res <- size.test.lc.mean.ws(.05, .90, 50.7, 2, .8, q)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 29)
})


test_that("size.equiv.mean.ps returns valid number", {
  res <- size.equiv.mean.ps(.10, .85, 15, .5, .7, 1.5)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 68)
})



test_that("size.supinf.mean.ps returns valid number", {
  res <- size.supinf.mean.ps(.05, .80, 225, 9, .75, 4)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 38)
})


test_that("size.test.mann returns valid number", {
  res <- size.test.mann(.05, .90, .3)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 44)
})


test_that("size.test.sign1 returns valid number", {
  res <- size.test.sign1(.05, .90, .3)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 56)
})


test_that("size.test.sign.ps returns valid number", {
  res <- size.test.sign.ps(.05, .90, .75)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 32)
})


test_that("size.test.cronbach returns valid number", {
  res <- size.test.cronbach(.05, .85, .80, 5, .7)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 139)
})


test_that("pi.score1 returns valid matrix", {
  colnames_expected <- c(
    "Predicted", "df", "LL", "UL"
  )
  
  res <- pi.score1(.05, 24.5, 3.65, 40)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("pi.score2 returns valid matrix", {
  colnames_expected <- c(
    "Predicted", "df", "LL", "UL"
  )
  
  res <- pi.score2(.05, 29.57, 18.35, 2.68, 1.92, 40, 45)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(2, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("random.sample returns valid vector", {

  res <- random.sample(3000, 25)
  
  testthat::expect_equal(length(res), 25)

})


test_that("randomize returns valid vector", {
  
  n <- c(10, 10, 5)
  res <- randomize(n)
  
  testthat::expect_equal(length(res), 25)
  
})


test_that("random.y returns valid vector", {
  
  n <- c(10, 10, 5)
  res <- random.y(10, 3.6, 2.8, 1, 7, 0) 
  
  testthat::expect_equal(length(res), 10)
  
})


test_that("ci.var.upper returns valid number", {
  res <- ci.var.upper(.25, 15, 60)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 17.2326356)
})


test_that("etasqr.adj returns valid number", {
  res <- etasqr.adj(.315, 2, 42)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 0.28238095)
})


test_that("pi.score2 returns valid matrix", {
  colnames_expected <- c(
    "F", "dfA", "dfE", "p", "eta-squared", "adj eta-squared"
  )
  
  m <- c(12.4, 8.6, 10.5)
  sd <- c(3.84, 3.12, 3.48)
  n <- c(20, 20, 20)
  res <- test.anova1.bs(m, sd, n)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})