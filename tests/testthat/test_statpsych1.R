library(statpsych)

test_that("ci.mean returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "LL", "UL")
  
  res <- ci.mean(.05, 24.5, 3.65, 40)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.stdmean returns valid matrix", {
  colnames_expected <- c("Estimate", "adj Estimate", "SE", "LL", "UL")
  
  res <- ci.stdmean(.05, 24.5, 3.65, 40, 20)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.mean2 returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "t", "df", "p", "LL", "UL")
  
  res <- ci.mean2(.05, 15.4, 10.3, 2.67, 2.15, 30, 20)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(2, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
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
  
  testthat::expect_snapshot(res)
})


test_that("ci.tukey returns valid matrix", {
  colnames_expected <- c("", "", "Estimate", "SE", "t", "df", "p", "LL", "UL")

  m <- c(12.86, 17.57, 26.29)
  sd <- c(3.185, 3.995, 3.773)
  n <- c(20, 20, 20)
  res <- ci.tukey(.05, m, sd, n)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(3, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.ratio.mean2  returns valid matrix", {
  colnames_expected <- c("Mean1", "Mean2", "Mean1/Mean2", "LL", "UL")
  
  y2 <- c(32, 39, 26, 35, 43, 27, 40, 37, 34, 29, 49, 42, 40)
  y1 <- c(36, 44, 47, 42, 49, 39, 46, 31, 33, 48)
  res <- ci.ratio.mean2(.05, y1, y2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.stdmean2  returns valid matrix", {
  colnames_expected <- c("Estimate", "adj Estimate", "SE", "LL", "UL")
  
  res <- ci.stdmean2(.05, 35.1, 26.7, 7.32, 6.98, 30, 30)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(4, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.stdmean.strat returns valid matrix", {
  colnames_expected <- c("Estimate", "adj Estimate", "SE", "LL", "UL")
  
  res <- ci.stdmean.strat(.05, 30.2, 30.8, 10.5, 11.2, 200, 200, .533)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(3, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.lc.stdmean.bs returns valid matrix", {
  colnames_expected <- c("Estimate", "adj Estimate", "SE", "LL", "UL")
  
  m <- c(33.5, 37.9, 38.0, 44.1)
  sd <- c(3.84, 3.84, 3.65, 4.98)
  n <- c(10,10,10,10)
  v <- c(.5, .5, -.5, -.5)
  res <- ci.lc.stdmean.bs(.05, m, sd, n, v)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(3, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.mean.ps returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "t", "df", "p", "LL", "UL")
  
  res <- ci.mean.ps(.05, 58.2, 51.4, 7.43, 8.92, .537, 30)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.ratio.mean.ps returns valid matrix", {
  colnames_expected <- c("Mean1", "Mean2", "Mean1/Mean2", "LL", "UL")
  
  y1 <- c(3.3, 3.6, 3.0, 3.1, 3.9, 4.2, 3.5, 3.3)
  y2 <- c(3.0, 3.1, 2.7, 2.6, 3.2, 3.8, 3.2, 3.0)
  res <- ci.ratio.mean.ps(.05, y1, y2)

  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.stdmean.ps returns valid matrix", {
  colnames_expected <- c("Estimate", "adj Estimate", "SE", "LL", "UL")
  
  res <- ci.stdmean.ps(.05, 110.4, 102.1, 15.3, 14.6, .75, 25)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(3, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.lc.stdmean.ws returns valid matrix", {
  colnames_expected <- c("Estimate", "adj Estimate", "SE", "LL", "UL")
  
  m <- c(33.5, 37.9, 38.0, 44.1)
  sd <- c(3.84, 3.84, 3.65, 4.98)
  q <- c(.5, .5, -.5, -.5)
  res <- ci.lc.stdmean.ws(.05, m, sd, .672, 20, q)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(2, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})



test_that("ci.mad returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "LL", "UL")
  
  y <- c(30, 20, 15, 10, 10, 60, 20, 25, 20, 30, 10, 5, 50, 40, 
         20, 10, 0, 20, 50)
  res <- ci.mad(.05, y)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.ratio.mad2 returns valid matrix", {
  colnames_expected <- c("MAD1", "MAD2", "MAD1/MAD2", "LL", "UL")
  
  y1 <- c(32, 39, 26, 35, 43, 27, 40, 37, 34, 29)
  y2 <- c(36, 44, 47, 42, 49, 39, 46, 31, 33, 48)
  res <- ci.ratio.mad2(.05, y1, y2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.ratio.mad.ps returns valid matrix", {
  colnames_expected <- c("MAD1", "MAD2", "MAD1/MAD2", "LL", "UL")
  
  y2 <- c(21, 4, 9, 12, 35, 18, 10, 22, 24, 1, 6, 8, 13, 16, 19)
  y1 <- c(67, 28, 30, 28, 52, 40, 25, 37, 44, 10, 14, 20, 28, 40, 51)
  res <- ci.ratio.mad.ps(.05, y1, y2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.cv returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "LL", "UL")
  
  res <- ci.cv(.05, 24.5, 3.65, 40)

  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.ratio.cv2 returns valid matrix", {
  colnames_expected <- c("Estimate", "LL", "UL")
  
  res <- ci.ratio.cv2(.05, 34.5, 26.1, 4.15, 2.26, 50, 50)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.cod returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "LL", "UL")
  
  y <- c(30, 20, 15, 10, 10, 60, 20, 25, 20, 30, 10, 5, 50, 40,
         20, 10, 0, 20, 50)
  res <- ci.cod(.05, y)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.median returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "LL", "UL")
  
  y <- c(30, 20, 15, 10, 10, 60, 20, 25, 20, 30, 10, 5, 50, 40,
         20, 10, 0, 20, 50)
  res <- ci.median(.05, y)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.median2 returns valid matrix", {
  colnames_expected <- c("Median1", "Median2", "Median1-Median2", "SE", "LL", "UL")
  
  y1 <- c(32, 39, 26, 35, 43, 27, 40, 37, 34, 29)
  y2 <- c(36, 44, 47, 42, 49, 39, 46, 31, 33, 48)
  res <- ci.median2(.05, y1, y2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.ratio.median2 returns valid matrix", {
  colnames_expected <- c("Median1", "Median2", "Median1/Median2", "LL", "UL")
  
  y2 <- c(32, 39, 26, 35, 43, 27, 40, 37, 34, 29, 49, 42, 40)
  y1 <- c(36, 44, 47, 42, 49, 39, 46, 31, 33, 48)
  res <- ci.ratio.median2(.05, y1, y2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
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
  
  testthat::expect_snapshot(res)
})


test_that("ci.median.ps returns valid matrix", {
  colnames_expected <- c(
    "Median1", "Median2", "Median1-Median2", "SE", "LL", "UL", "SE1", "SE2", "COV"
  )
  
  y1 <- c(21, 4, 9, 12, 35, 18, 10, 22, 24, 1, 6, 8, 13, 16, 19)
  y2 <- c(67, 28, 30, 28, 52, 40, 25, 37, 44, 10, 14, 20, 28, 40, 51)
  res <- ci.median.ps(.05, y1, y2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
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
  
  testthat::expect_snapshot(res)
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
  
  testthat::expect_snapshot(res)
})


test_that("ci.random.anova returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "LL", "UL"
  )
  
  m <- c(56.1, 51.2, 60.3, 68.2, 48.9, 70.5)
  sd <- c(9.45, 8.79, 9.71, 8.90, 8.31, 9.75)
  res <- ci.random.anova(.05, m, sd, 20)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(4, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.cronbach returns valid matrix", {
  colnames_expected <- c(
   "Estimate", "SE", "LL", "UL"
  )
  
  res <- ci.cronbach(.05, .85, 7, 89)
  
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("size.ci.mean returns valid number", {
  colnames_expected <- c(
    "Sample size"
  )
  
  
  res <- size.ci.mean(.05, 264.4, 10)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  testthat::expect_equal(res[[1,1]], 43)
  
  testthat::expect_snapshot(res)
})


test_that("size.ci.mean2 returns valid matrix", {
  colnames_expected <- c(
    "n1", "n2"
  )
  
  res <- size.ci.mean2(.05, 37.1, 5, 1)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("size.ci.stdmean2 returns valid matrix", {
  colnames_expected <- c(
    "n1", "n2"
  )
  
  res <- size.ci.stdmean2(.05, .75, .5, 1)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(2, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("size.ci.ratio.mean2 returns valid matrix", {
  colnames_expected <- c(
    "n1", "n2"
  )
  
  res <- size.ci.ratio.mean2(.05, .4, 3.5, 3.1, 1.2, 2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("size.ci.lc.mean.bs returns valid number", {

  v <- c(.5, .5, -1)
  res <- size.ci.lc.mean.bs(.05, 5.62, 2.0, v)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1]], 34)
  
  testthat::expect_snapshot(res)
})


test_that("size.ci.stdmean.ps returns valid number", {
  res <- size.ci.mean.ps(.05, 265, .8, 10)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1]], 19)
  
  testthat::expect_snapshot(res)
})


test_that("size.ci.ratio.mean2 returns valid number", {
  res <- size.ci.stdmean.ps(.05, 1, .65, .6)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1]], 46)
  
  testthat::expect_snapshot(res)
})


test_that("size.ci.ratio.mean.ps returns valid number", {
  res <- size.ci.ratio.mean.ps(.05, 400, 150, 100, .7, 1.2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1]], 21)
  
  testthat::expect_snapshot(res)
})


test_that("size.ci.lc.stdmean.ws returns valid matrix", {
  q <- c(.5, .5, -.5, -.5)
  res <- size.ci.lc.mean.ws(.05, 265, .8, 10, q)

  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 11)
  
  testthat::expect_snapshot(res)
})


test_that("size.ci.lc.mean.ws returns valid matrix", {
  q <- c(.5, .5, -.5, -.5)
  res <- size.ci.lc.stdmean.ws(.05, 1, .7, .6, q)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 26)
  
  testthat::expect_snapshot(res)
})


test_that("size.ci.cronbach returns valid number", {
  res <- size.ci.cronbach(.05, .85, 5, .1)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 89)
  
  testthat::expect_snapshot(res)
})


test_that("size.ci.second returns valid number", {
  res <- size.ci.second(20, 5.3, 2.5)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 70)
  
  testthat::expect_snapshot(res)
})


test_that("size.test.mean returns valid number", {
  res <- size.test.mean(.05, .9, 80.5, 7)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 20)
  
  testthat::expect_snapshot(res)
})


test_that("size.test.mean2 returns valid matrix", {
  colnames_expected <- c(
    "n1", "n2"
  )
  
  res <- size.test.mean2(.05, .95, 100, 10, 1) 
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("size.test.lc.mean.bs returns valid matrix", {
  v <- c(1, -1, -1, 1)
  res <- size.test.lc.mean.bs(.05, .90, 27.5, 5, v)

  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 47)  
  
  testthat::expect_snapshot(res)
    
})


test_that("size.equiv.mean2 returns valid matrix", {

  res <- size.equiv.mean2(.10, .80, 15, 2, 4)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 50)
  
  testthat::expect_snapshot(res)
})


test_that("size.supinf.mean2 returns valid matrix", {
  res <- size.supinf.mean2(.05, .80, 225, 9, 4)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 143)
  
  testthat::expect_snapshot(res)
})


test_that("size.test.mean.ps returns valid number", {
  res <- size.test.mean.ps(.05, .80, 1.25, .5, .75) 
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 22)
  
  testthat::expect_snapshot(res)
})


test_that("size.test.lc.mean.ws returns valid matrix", {
  q <- c(.5, .5, -.5, -.5)
  res <- size.test.lc.mean.ws(.05, .90, 50.7, 2, .8, q)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 29)
  
  testthat::expect_snapshot(res)
})


test_that("size.equiv.mean.ps returns valid number", {
  res <- size.equiv.mean.ps(.10, .85, 15, .5, .7, 1.5)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 68)
  
  testthat::expect_snapshot(res)
})



test_that("size.supinf.mean.ps returns valid number", {
  res <- size.supinf.mean.ps(.05, .80, 225, 9, .75, 4)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 38)
  
  testthat::expect_snapshot(res)
})


test_that("size.test.mann returns valid number", {
  res <- size.test.mann(.05, .90, .3)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 44)
  
  testthat::expect_snapshot(res)
})


test_that("size.test.sign returns valid number", {
  res <- size.test.sign(.05, .90, .3)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 67)
  
  testthat::expect_snapshot(res)
})


test_that("size.test.sign.ps returns valid number", {
  res <- size.test.sign.ps(.05, .90, .75)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 42)
  
  testthat::expect_snapshot(res)
})


test_that("size.test.cronbach returns valid number", {
  res <- size.test.cronbach(.05, .85, .80, 5, .7)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 139)
  
  testthat::expect_snapshot(res)
})


test_that("pi.score returns valid matrix", {
  colnames_expected <- c(
    "Predicted", "df", "LL", "UL"
  )
  
  res <- pi.score(.05, 24.5, 3.65, 40)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("pi.score2 returns valid matrix", {
  colnames_expected <- c(
    "Predicted", "df", "LL", "UL"
  )
  
  res <- pi.score2(.05, 29.57, 18.35, 2.68, 1.92, 40, 45)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(2, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
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
  
  testthat::expect_snapshot(res)
})


test_that("etasqr.adj returns valid number", {
  res <- etasqr.adj(.315, 2, 42)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 0.28238095)
  
  testthat::expect_snapshot(res)
})


test_that("test.anova.bs returns valid matrix", {
  colnames_expected <- c(
    "F", "dfA", "dfE", "p", "Eta-squared", "adj Eta-squared"
  )
  
  m <- c(12.4, 8.6, 10.5)
  sd <- c(3.84, 3.12, 3.48)
  n <- c(20, 20, 20)
  res <- test.anova.bs(m, sd, n)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("etasqr.gen.2way returns valid matrix", {
  colnames_expected <- c(
    "A", "B", "AB"
  )
  
  res <- etasqr.gen.2way(12.3, 15.6, 5.2, 7.9)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(3, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.ratio.cod2 returns valid matrix", {
  colnames_expected <- c(
    "COD1",      "COD2", "COD1/COD2",       "LL",       "UL"
  )

  y1 <- c(32, 39, 26, 35, 43, 27, 40, 37, 34, 29)
  y2 <- c(36, 44, 47, 42, 49, 39, 46, 31, 33, 48)
  res <- ci.ratio.cod2(.05, y1, y2)

  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.etasqr returns valid matrix", {
  colnames_expected <- c(
    "Eta-squared", "adj Eta-squared",     "SE",   "LL",        "UL"
  )
  
  res <- ci.etasqr(.05, .241, 3, 116)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.reliability returns valid vector", {
  colnames_expected <- c("Estimate", "LL", "UL")
  
  
  res <- ci.reliability(.05, .88, .0147, 100)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, 3))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)

})


test_that("ci.sign returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE",        "LL",        "UL"
  )
  
  y <- c(30, 20, 15, 10, 10, 60, 20, 25, 20, 30, 10, 5, 50, 40, 20, 10, 0, 20, 50)
  res <- ci.sign(.05, y, 9)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.slope.mean.bs returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "t", "df", "p",       "LL",        "UL"
  )
  
  m <- c(33.5, 37.9, 38.0, 44.1)
  sd <- c(3.84, 3.84, 3.65, 4.98)
  n <- c(10,10,10,10)
  x <- c(5, 10, 20, 30)
  res <- ci.slope.mean.bs(.05, m, sd, n, x)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(2, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.etasqr returns valid matrix", {
  colnames_expected <- c(
    "Kurtosis", "Excess", "p"
  )
  
  y <- c(30, 20, 15, 10, 10, 60, 20, 25, 20, 30, 10, 5, 50, 40, 95)
  res <- test.kurtosis(y)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  # testthat::expect_snapshot(res)
})



test_that("test.skew returns valid matrix", {
  colnames_expected <- c(
    "Skewness", "p"
  )
  
  y <- c(30, 20, 15, 10, 10, 60, 20, 25, 20, 30, 10, 5, 50, 40, 95)
  res <- test.skew(y)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  # testthat::expect_snapshot(res)
})


test_that("test.mono.mean.bs returns valid matrix", {
  colnames_expected <- c(
    "", "", "Estimate", "SE",        "LL",        "UL"
  )
  
  m <- c(12.86, 24.57, 36.29, 53.21)
  sd <- c(13.185, 12.995, 14.773, 15.145)
  n <- c(20, 20, 20, 20)
  res <- test.mono.mean.bs(.05, m, sd, n)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(3, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})

# Simulation tests commented out because they take too long for CRAN

# test_that("sim.ci.median2 returns valid matrix", {
#   colnames_expected <- c(
#     "Coverage", "Lower Error", "Upper Error", "Ave CI Width"
#   )
#   
#   res <- sim.ci.median2(.05, 20, 20, 2, 5, 4, 5000)
#   
#   testthat::expect_equal(class(res), c("matrix", "array"))
#   testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
#   testthat::expect_equal(colnames(res), colnames_expected)
# })

# test_that("sim.ci.median.ps returns valid matrix", {
#   colnames_expected <- c(
#     "Coverage", "Lower Error", "Upper Error", "Ave CI Width"
#   )
#   
#   res <- sim.ci.median.ps(.05, 30, 1.5, .7, 4, 3, 2000)
#   
#   testthat::expect_equal(class(res), c("matrix", "array"))
#   testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
#   testthat::expect_equal(colnames(res), colnames_expected)
# })
# 
# test_that("sim.ci.stdmean2 returns valid matrix", {
#   colnames_expected <- c(
#     "Coverage", "Lower Error", "Upper Error", "Ave CI Width", "Ave Est"
#   )
#   
#   res <- sim.ci.stdmean2(.05, 20, 20, 1.5, 3, 4, .75, 5000)
#   
#   testthat::expect_equal(class(res), c("matrix", "array"))
#   testthat::expect_equal(dim(res), c(2, length(colnames_expected)))
#   testthat::expect_equal(colnames(res), colnames_expected)
# })

test_that("power.mean returns valid matrix", {
  colnames_expected <- c(
    "Power"
  )
  
  res <- power.mean(.05, 15, 80.5, 7)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("power.mean2 returns valid matrix", {
  colnames_expected <- c(
    "Power"
  )
  
  res <- power.mean2(.05, 25, 25, 5.0, 6.0, 2)

  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("power.mean.ps returns valid matrix", {
  colnames_expected <- c(
    "Power"
  )
  
  res <- power.mean.ps(.05, 20, 10.0, 12.0, 2, .7)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("power.lc.bs returns valid matrix", {
  colnames_expected <- c(
    "Power"
  )

  n <- c(20, 20, 20, 20)
  var <- c(70, 70, 80, 80)
  v <- c(.5, .5, -.5, -.5)
  res <- power.lc.mean.bs(.05, n, var, 5, v)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})



test_that("ci.cqv returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE",        "LL",        "UL"
  )
  
  y <- c(30, 20, 15, 10, 10, 60, 20, 25, 20, 30, 10, 5, 50, 40,
         20, 10, 0, 20, 50)
  res <- ci.cqv(.05, y)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.ratio.sd2 returns valid matrix", {
  colnames_expected <- c(
    "SD1", "SD2", "SD1/SD2", "LL",        "UL"
  )
  
  y1 <- c(32, 39, 26, 35, 43, 27, 40, 37, 34, 29)
  y2 <- c(36, 44, 47, 42, 49, 39, 46, 31, 33, 48)
  res <- ci.ratio.sd2(.05, y1, y2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("size.ci.etasqr returns valid matrix", {
  colnames_expected <- c(
    "Sample size per group"
  )
  
  res <- size.ci.etasqr(.05, .333, 3, .2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.2x2.stdmean.bs returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "adj Estimate", "SE", "LL",        "UL"
  )
  
  y11 = c(14, 15, 11, 7, 16, 12, 15, 16, 10, 9)
  y12 = c(18, 24, 14, 18, 22, 21, 16, 17, 14, 13)
  y21 = c(16, 11, 10, 17, 13, 18, 12, 16, 6, 15)
  y22 = c(18, 17, 11, 9, 9, 13, 18, 15, 14, 11)
  res <- ci.2x2.stdmean.bs(.05, y11, y12, y21, y22)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(7, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.2x2.median.bs returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL",        "UL"
  )
  
  y11 <- c(14, 15, 11, 7, 16, 12, 15, 16, 10, 9)
  y12 <- c(18, 24, 14, 18, 22, 21, 16, 17, 14, 13)
  y21 <- c(16, 11, 10, 17, 13, 18, 12, 16, 6, 15)
  y22 <- c(18, 17, 11, 9, 9, 13, 18, 15, 14, 11)
  res <- ci.2x2.median.bs(.05, y11, y12, y21, y22)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(7, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.2x2.stdmean.ws returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "adj Estimate",  "SE", "LL",        "UL"
  )
  
  y11 <- c(21, 39, 32, 29, 27, 17, 27, 21, 28, 17, 12, 27)
  y12 <- c(20, 36, 33, 27, 28, 14, 30, 20, 27, 15, 11, 22)
  y21 <- c(21, 36, 30, 27, 28, 15, 27, 18, 29, 16, 11, 22)
  y22 <- c(18, 34, 29, 28, 28, 17, 27, 21, 26, 16, 14, 23)
  res <- ci.2x2.stdmean.ws(.05, y11, y12, y21, y22)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(7, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})



test_that("ci.2x2.stdmean.mixed returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "adj Estimate",  "SE", "LL",        "UL"
  )
  
  y11 <- c(18, 19, 20, 17, 20, 16)
  y12 <- c(19, 18, 19, 20, 17, 16)
  y21 <- c(19, 16, 16, 14, 16, 18)
  y22 <- c(16, 10, 12,  9, 13, 15)
  res <- ci.2x2.stdmean.mixed(.05, y11, y12, y21, y22)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(7, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
  
})


test_that("ci.2x2.median.mixed returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL",        "UL"
  )
  
  y11 <- c(18, 19, 20, 17, 20, 16)
  y12 <- c(19, 18, 19, 20, 17, 16)
  y21 <- c(19, 16, 16, 14, 16, 18)
  y22 <- c(16, 10, 12,  9, 13, 15)
  res <- ci.2x2.median.mixed(.05, y11, y12, y21, y22)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(7, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.2x2.median.w returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL",        "UL"
  )
  
  y11 <- c(221, 402, 333, 301, 284, 182, 281, 230, 290, 182, 133, 278)
  y12 <- c(221, 371, 340, 288, 293, 150, 317, 211, 286, 161, 126, 234)
  y21 <- c(219, 371, 314, 279, 284, 155, 278, 185, 296, 169, 118, 229)
  y22 <- c(170, 332, 280, 273, 272, 160, 260, 204, 252, 153, 137, 221)
  res <- ci.2x2.median.ws(.05, y11, y12, y21, y22)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(7, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("pi.var example", {
  res <- pi.var(.05, 15, 40, 100, 2)
  testthat::expect_snapshot(res)
})


test_that("ci.bayes.normal returns valid matrix", {
  colnames_expected <- c(
    "Posterior mean", "Posterior SD", "LL",        "UL"
  )
  
  res <- ci.bayes.normal(.05, 30, 2, 24.5, 0.577)

  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("spearmanbrown returns valid matrix", {
  colnames_expected <- c(
    "Reliability of r2 measurements"
  )
  
  res <- spearmanbrown(.6, 10, 20)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})



test_that("ci.mean.fpc returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL",        "UL"
  )
  
  res <- ci.mean.fpc(.05, 24.5, 3.65, 40, 300)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("test.mean returns valid number", {
  colnames_expected <- c(
    "t", "df", "p"
  )
  
  
  res <- test.mean(24.5, 3.65, 40, 23)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})



test_that("size.ci.mean.prior returns valid matrix", {
  colnames_expected <- c(
    "Sample size"
  )
  
  res <- size.ci.mean.prior(.05, .10, 26.4, 25, 4)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(res[[1,1]], 44)
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})

