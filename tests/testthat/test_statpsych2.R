library(statpsych)

test_that("ci.cor returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  res <- ci.cor(.05, .536, 0, 50)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.spcor returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  res <- ci.spcor(.05, .582, .699, 20)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.cor2 returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  res <- ci.cor2(.05, .886, .802, 200, 200)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.cor.dep returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  res <- ci.cor.dep(.05, .396, .179, .088, 166)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.cor2.gen returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "LL", "UL"
  )
  
  res <- ci.cor2.gen(.4, .35, .47, .2, .1, .32)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.pbcor returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  res <- ci.pbcor(.05, 28.32, 21.48, 3.81, 3.09, 40, 40)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(2, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.spear returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  y <- c(21, 4, 9, 12, 35, 18, 10, 22, 24, 1, 6, 8, 13, 16, 19)
  x <- c(67, 28, 30, 28, 52, 40, 25, 37, 44, 10, 14, 20, 28, 40, 51)
  res <- ci.spear(.05, y, x)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.spear2 returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  y <- c(21, 4, 9, 12, 35, 18, 10, 22, 24, 1, 6, 8, 13, 16, 19)
  x <- c(67, 28, 30, 28, 52, 40, 25, 37, 44, 10, 14, 20, 28, 40, 51)
  res <- ci.spear2(.05, .54, .48, 180, 200)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.mape returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  res <- c(-2.70, -2.69, -1.32, 1.02, 1.23, -1.46, 2.21, -2.10, 2.56,
         -3.02, -1.55, 1.46, 4.02, 2.34)
  res <- ci.mape(.05, res, 1)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.condslope returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "t", "df", "p", "LL", "UL"
  )
  
  res <- ci.condslope(.05, .132, .154, .031, .021, .015, 5.2, 10.6, 122)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(2, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.lc.reg returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "t", "df", "p", "LL", "UL"
  )
  
  est <- c(1.74, 1.83, 0.482)
  se <- c(.483, .421, .395)
  n <- c(40, 40, 40)
  v <- c(.5, .5, -1)
  res <- ci.lc.reg(.05, est, se, n, 4, v)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.fisher returns valid matrix", {
  colnames_expected <- c(
   "Estimate", "LL", "UL"
  )
  
  res <- ci.fisher(.05, .641, .052)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.indirect returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  res <- ci.indirect (.05, 2.48, 1.92, .586, .379)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  #testthat::expect_snapshot(res)
})


test_that("size.ci.slope returns valid numeric", {
  
  x <- c(2, 5, 8)
  res <- size.ci.slope(.05, 31.1, x, 1)

  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 83)
  
  testthat::expect_snapshot(res)
})


test_that("size.ci.cor returns valid numeric", {
  
  res <- size.ci.cor(.05, .362, 0, .25)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 188)
  
  testthat::expect_snapshot(res)
})


test_that("size.ci.spear returns valid numeric", {
  
  colnames_expected <- c(
    "Sample size"
  )
  
  res <- size.ci.spear(.05, .362, .25)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 200)
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("size.ci.pbcor returns valid numeric", {
  colnames_expected <- c(
    "Sample size"
  )

  res <- size.ci.pbcor(.05, .40, .25, .73)

  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 168)
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("size.ci.pbcor returns valid numeric", {
  colnames_expected <- c(
    "Sample size"
  )
  
  res <- size.ci.pbcor(.05, .40, .25, .73)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 168)
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("size.ci.rsqr returns valid numeric", {
  
  res <- size.ci.rsqr(.05, .333, 2, .2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 226)
  
  testthat::expect_snapshot(res)
})


test_that("size.ci.condmean returns valid numeric", {
  
  res <- size.ci.condmean(.05, 120, 125, 15, 5)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 210)
  
  testthat::expect_snapshot(res)
})


test_that("size.ci.lc.ancova returns valid numeric", {
  
  v <- c(1, -1)
  res <- size.ci.lc.ancova(.05, 1.37, 1, 0, 1.5, v)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 21)
  
  testthat::expect_snapshot(res)
})


test_that("size.test.slope returns valid numeric", {
  
  x <- c(2, 5, 8)
  res <- size.test.slope(.05, .9, 31.1, x, .75, 0)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 100)
  
  testthat::expect_snapshot(res)
})


test_that("size.test.cor returns valid numeric", {
  
  res <- size.test.cor(.05, .9, .45, 0, 0)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 48)
  
  testthat::expect_snapshot(res)
})



test_that("size.interval.cor returns valid numeric", {
  
  res <- size.interval.cor(.05, .8, .1, 0, .25)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 360)
  
  testthat::expect_snapshot(res)
})


test_that("size.test.lc.ancova returns valid numeric", {
  
  v <- c(.5, .5, -1)
  res <- size.test.lc.ancova(.05, .9, 1.37, .7, 1, 0, v)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 46)
  
  testthat::expect_snapshot(res)
})


test_that("ci.indirect returns valid matrix", {
  colnames_expected <- c("Coefficient")
  
  x <- c(25, 50, 75, 100)
  res <- slope.contrast(x)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(4, 1))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  # testthat::expect_snapshot(res)
})
  

test_that("random.yx returns valid data.frame", {
  colnames_expected <- c(
    "y", "x"
  )
  
  res <- random.yx(10, 50, 20, 4, 2, .5, 1)
  
  testthat::expect_equal(class(res), c("data.frame"))
  testthat::expect_equal(dim(res), c(10, length(res)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  #testthat::expect_snapshot(res)
})


test_that("ci.rsqr returns valid matrix", {
  colnames_expected <- c(
    "R-squared",    "adj R-squared",  "SE",        "LL",        "UL"
  )
  
  res <- ci.rsqr(.05, .241, 3, 116)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.lc.gen.bs returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  est <- c(3.86, 4.57, 2.29, 2.88)
  se <- c(0.185, 0.365, 0.275, 0.148)
  v <- c(.5, .5, -.5, -.5)
  res <- ci.lc.gen.bs(.05, est, se, v)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.lc.glm returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "t", "df", "p", "LL", "UL"
  )
  
  y <- c(43, 62, 49, 60, 36, 79, 55, 42, 67, 50)
  x1 <- c(3, 6, 4, 6, 2, 7, 4, 2, 7, 5)
  x2 <- c(4, 6, 3, 7, 1, 9, 3, 3, 8, 4)
  out <- lm(y ~ x1 + x2)
  b <- coef(out)
  V <- vcov(out)
  n <- length(y)
  q <- c(0, .5, .5)
  b
  res <- ci.lc.glm(.05, n, b, V, q)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.theil returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  y <- c(21, 4, 9, 12, 35, 18, 10, 22, 24, 1, 6, 8, 13, 16, 19)
  x <- c(67, 28, 30, 28, 52, 40, 25, 37, 44, 10, 14, 20, 28, 40, 51)
  res <- ci.theil(.05, y, x)

  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("power.cor returns valid matrix", {
  colnames_expected <- c(
    "Power"
  )
  
  res <- power.cor(.05, 80, .3, 0, 0)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})

test_that("power.cor2 returns valid matrix", {
  colnames_expected <- c(
    "Power"
  )
  
  res <- power.cor2(.05, 200, 200, .4, .2, 0)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("size.test.cor2 returns valid matrix", {
  colnames_expected <- c(
    "n1", "n2"
  )
  
  res <- size.test.cor2(.05, .8, .4, .2, 0, 1)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  testthat::expect_equal(res[[1,1]], 325)
  testthat::expect_equal(res[[1,2]], 325)
  
  res <- size.test.cor2(.05, .8, .4, .2, 0, 2)
  testthat::expect_equal(res[[1,1]], 245)
  testthat::expect_equal(res[[1,2]], 490)
  
  testthat::expect_snapshot(res)
  
})


test_that("size.test.cronbach2 returns valid matrix", {
  colnames_expected <- c(
    "Sample size per group"
  )
  
  res <- size.test.cronbach2(.05, .80, .85, .70, 8)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("size.ci.cronbach2 returns valid matrix", {
  colnames_expected <- c(
    "Sample size per group"
  )
  
  res <- size.ci.cronbach2(.05, .85, .70, 8, .15)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("size.ci.mape returns valid numeric", {
  colnames_expected <- c(
    "Sample size"
  )
  
  res <- size.ci.mape(.05, 4.5, 5, 2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 57)
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("size.ci.indirect returns valid matrix", {
  colnames_expected <- c(
    "Sample size"
  )
  
  res <- size.ci.indirect(.05, .4, .5, .2)

  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.ratio.mape2 returns valid matrix", {
  colnames_expected <- c(
    "MAPE1", "MAPE2", "MAPE1/MAPE2", "LL", "UL"
  )

  res1 <- c(-2.70, -2.69, -1.32, 1.02, 1.23, -1.46, 2.21, -2.10, 2.56, -3.02
          -1.55, 1.46, 4.02, 2.34)
  res2 <- c(-0.71, -0.89, 0.72, -0.35, 0.33 -0.92, 2.37, 0.51, 0.68, -0.85,
          -0.15, 0.77, -1.52, 0.89, -0.29, -0.23, -0.94, 0.93, -0.31 -0.04)
  res <- ci.ratio.mape2(.05, res1, res2, 1, 1)

  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.rel2 returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "LL", "UL"
  )
  
  res <- ci.rel2(.4, .35, .47, .2, .1, .32)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.cronbach2 returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "LL", "UL"
  )
  
  res <- ci.cronbach2(.05, .88, .76, 8, 8, 200, 250)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("ci.bscor returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  res <- ci.bscor(.05, 28.32, 21.48, 3.81, 3.09, 40, 40)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("pi.cor returns valid matrix", {
  colnames_expected <- c(
    "LL", "UL"
  )
  
  res <- pi.cor(.1, .761, 50, 100, 1)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  testthat::expect_equal(res[1,1], 0.60340924)
  testthat::expect_equal(res[1,2], 0.85732237)
  
  testthat::expect_snapshot(res)
  
})



test_that("test.cor returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "z", "p"
  )
  
  res <- test.cor(.484, 100, 0, .2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("test.spear returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "z", "p"
  )
  
  res <- test.spear(.471, .2, 100)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("test.cor2 returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "z", "p"
  )
  
  res <- test.cor2(.684, .437, 100, 125, 0)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("test.spear2 returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "z", "p"
  )
  
  res <- test.spear2(.684, .437, 100, 125)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("size.ci.cor2 returns valid matrix", {
  colnames_expected <- c(
    "Sample size per group"
  )
  
  res <- size.ci.cor2(.05, .8, .5, .2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(res[[1,1]], 271)
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("size.ci.spear2 returns valid matrix", {
  colnames_expected <- c(
    "Sample size per group"
  )
  
  res <- size.ci.spear2(.05, .8, .5, .2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(res[[1,1]], 314)
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("size.ci.cor.prior returns valid matrix", {
  colnames_expected <- c(
    "Sample size"
  )
  
  res <- size.ci.cor.prior(.05, .10, .438, 100, .2)
         
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(res[[1,1]], 331)
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("adj.se returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "adj SE", "t", "df", "p", "LL", "UL"
  )
  
  se <- c(1.57, 3.15, 0.982)
  b <- c(3.78, 8.21, 2.99)
  res <- adj.se(.05, 10.26, 8.37, 114, se, b)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(3, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("fitindices returns valid matrix", {
  colnames_expected <- c(
    "NFI", "adj NFI", "CFI", "TLI", "RMSEA"
  )
  
  res <- fitindices(14.21, 10, 258.43, 20, 300)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})



test_that("size.ci.biphi example", {
  res <- size.ci.biphi(.05, .2, .5, .3, .4)
  testthat::expect_snapshot(res)
})



test_that("size.ci.ancova2 example", {
  res <- size.ci.ancova2(.05, 1.37, 1, 0, 1.5, 2)
  testthat::expect_snapshot(res)
})


test_that("size.ci.slope.gen example", {
  res <- size.ci.slope.gen(.05, 3.15, 50, 5)
  testthat::expect_snapshot(res)
})


test_that("size.test.ancova2 example", {
  res <- size.test.ancova2(.05, .9, 1.37, .7, 1, 0, 2)
  testthat::expect_snapshot(res)
})


test_that("size.test.slope.gen example", {
  res <- size.test.slope.gen(.05, .8, 3.15, 50, 5)
  testthat::expect_snapshot(res)
})


test_that("signal example", {
  res <- signal(82, 46, 100, 100)
  testthat::expect_snapshot(res)
})


test_that("expon.slope example", {
  res <- expon.slope(.05, .502, .0396)
  testthat::expect_snapshot(res)
})



test_that("ci.bayes.cor example", {
  res <- ci.bayes.cor(.05, .1, .536, 0, 50)
  testthat::expect_snapshot(res)
})


test_that("ci.bayes.spcor example", {
  res <- ci.bayes.spcor(.05, .1, .582, .137)
  testthat::expect_snapshot(res)
})


test_that("test.mono.median.bs example", {
  m <- c(12.86, 24.57, 36.29, 53.21)
  se <- c(2.85, 2.99, 3.73, 3.88)
  res <- test.mono.median.bs(.05, m, se)
  testthat::expect_snapshot(res)
})


test_that("ci.slope.median.bs example", {
  m <- c(33.5, 37.9, 38.0, 44.1)
  se <- c(0.84, 0.94, 1.65, 2.98)
  x <- c(5, 10, 20, 30)
  res <- ci.slope.median.bs(.05, m, se, x)
  testthat::expect_snapshot(res)
})


test_that("size.ci.gen example", {
  res <- size.ci.gen(.05, 2.89, 30, 8)
  testthat::expect_snapshot(res)
})


test_that("size.ci.gen2 example", {
  res <- size.ci.gen2(.05, .175, 30, .8, 1)
  testthat::expect_snapshot(res)
})


test_that("size.test.gen(.05, .8, 2.89, 30, 5) example", {
  res <- size.test.gen(.05, .8, 2.89, 30, 5)
  testthat::expect_snapshot(res)
})


test_that("size.test.gen2(.05, .85, .175, 30, .5, 1) example", {
  res <- size.test.gen2(.05, .85, .175, 30, .5, 1)
  testthat::expect_snapshot(res)
})