library(statpsych)


test_that("ci.prop returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  res <- ci.prop(.05, 12, 100)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(2, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


# test_that("ci.pairs.prop returns valid matrix", {
#   colnames_expected <- c(
#     "", "", "Estimate", "SE", "LL", "UL"
#   )
#   
#   f <- c(125, 82, 92)
#   res <- ci.pairs.prop(.05, f)
#   
#   testthat::expect_equal(class(res), c("matrix", "array"))
#   testthat::expect_equal(dim(res), c(3, length(colnames_expected)))
#   testthat::expect_equal(colnames(res), colnames_expected)
# })


test_that("ci.prop2 returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  res <- ci.prop2(.05, 35, 21, 150, 150)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.ratio.prop2 returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "LL", "UL"
  )
  
  res <- ci.ratio.prop2(.05, 35, 21, 150, 150)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.lc.prop.bs returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "z", "p", "LL", "UL"
  )
  
  f <- c(26, 24, 38)
  n <- c(60, 60, 60)
  c <- c(-.5, -.5, 1)
  res <- ci.lc.prop.bs(.05, f, n, c)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.pairs.prop.bs returns valid matrix", {
  colnames_expected <- c(
    "", "", "Estimate", "SE", "z", "p", "LL", "UL"
  )
  
  f <- c(111, 161, 132)
  n <- c(200, 200, 200)
  res <- ci.pairs.prop.bs(.05, f, n)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(3, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})



test_that("ci.slope.prop.bs returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "z", "p", "LL", "UL"
  )
  
  f <- c(14, 27, 38)
  n <- c(100, 100, 100)
  x <- c(10, 20, 40)
  res <- ci.slope.prop.bs(.05, f, n, x)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.prop.ps returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  res <- ci.prop.ps(.05, 12, 26, 4, 6)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.ratio.prop.ps returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "LL", "UL"
  )
  
  res <- ci.ratio.prop.ps(.05, 12, 26, 4, 6)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.condslope.log returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "exp(Estimate)", "z", "p", "LL", "UL"
  )
  
  res <- ci.condslope.log(.05, .132, .154, .031, .021, .015, 5.2, 10.6)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(2, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.oddsratio returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  res <- ci.oddsratio(.05, 229, 28, 96, 24)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.yule returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  res <- ci.yule(.05, 229, 28, 96, 24)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(4, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.phi returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  res <- ci.phi(.05, 229, 28, 96, 24)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.biphi returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  res <- ci.biphi(.05, 46, 15, 100, 100)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.tetra returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  res <- ci.tetra(.05, 46, 15, 54, 85)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.kappa returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  res <- ci.kappa(.05, 31, 12, 4, 58)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(2, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.agree returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  res <- ci.agree(.05, 100, 80, 4)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.popsize returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL", "UL"
  )
  
  res <- ci.popsize(.05, 794, 710, 741)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("test.prop returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "z", "p"
  )
  
  res <- test.prop(9, 20, .2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("test.prop2 returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "z", "p"
  )
  
  res <- test.prop2(11, 26, 50, 50)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("test.prop.bs returns valid matrix", {
  colnames_expected <- c(
    "Chi-square", "df", "p"
  )
  
  f <- c(35, 30, 15)
  n <- c(50, 50, 50)
  res <- test.prop.bs (f, n)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("test.prop.ps returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "z", "p"
  )
  
  res <- test.prop.ps(156, 96, 68, 80)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("size.ci.prop returns valid numeric", {

  res <- size.ci.prop(.05, .4, .2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 93)
})


test_that("size.ci.prop2 returns valid numeric", {
  
  res <- size.ci.prop2(.05, .4, .2, .15)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 274)
})


test_that("size.ci.ratio.prop2 returns valid numeric", {
  
  res <- size.ci.ratio.prop2(.05, .2, .1, 2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 416)
})


test_that("size.ci.lc.prop.bs returns valid numeric", {
  
  p <- c(.25, .30, .50, .50)
  v <- c(.5, .5, -.5, -.5)
  res <- size.ci.lc.prop.bs(.05, p, .2, v)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 87)
})


test_that("size.ci.prop.ps returns valid numeric", {
  
  p <- c(.25, .30, .50, .50)
  v <- c(.5, .5, -.5, -.5)
  res <- size.ci.prop.ps(.05, .2, .3, .8, .1)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 118)
})


test_that("size.ci.ratio.prop.ps returns valid numeric", {
  
  res <- size.ci.ratio.prop.ps(.05, .4, .2, .7, 2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 67)
})


test_that("size.ci.agree returns valid numeric", {
  
  res <- size.ci.agree(.05, .8, .2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 139)
})


test_that("size.test.prop returns valid numeric", {
  
  res <- size.test.prop(.05, .9, .5, .3)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 67)
})


test_that("size.test.prop2 returns valid numeric", {
  
  res <- size.test.prop2(.05, .8, .5, .5, .2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 109)
})


test_that("size.test.lc.prop.bs returns valid numeric", {
  
  p <- c(.25, .30, .50, .50)
  v <- c(.5, .5, -.5, -.5)
  res <- size.test.lc.prop.bs(.05, .9, p, .15, v)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 105)
})


test_that("size.equiv.prop2 returns valid numeric", {
  
  res <- size.equiv.prop2(.1, .8, .30, .35, .15)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 288)
})


test_that("size.supinf.prop2 returns valid numeric", {
  
  res <- size.supinf.prop2(.05, .9, .35, .20, .05)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 408)
})


test_that("size.test.prop.ps returns valid numeric", {
  
  res <- size.test.prop.ps(.05, .80, .4, .3, .5, .1)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 177)
})



test_that("size.equiv.prop.ps returns valid numeric", {
  
  res <- size.equiv.prop.ps(.1, .8, .30, .35, .40, .15)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 173)
})


test_that("size.supinf.prop.ps returns valid numeric", {
  
  res <- size.supinf.prop.ps(.05, .9, .35, .20, .45, .05)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(res[[1,1]], 227)
})


test_that("iqv returns valid matrix", {
  colnames_expected <- c(
    "Simpson", "Berger", "Shannon"
  )
  
  f <- c(10, 46, 15, 3)
  res <- iqv(f)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})



test_that("test.mono.prop.bs returns valid matrix", {
  colnames_expected <- c(
    "", "", "Estimate",         "SE",         "LL",        "UL"
  )
  
  f <- c(67, 49, 30, 10)
  n <- c(100, 100, 100, 100)
  res <- test.mono.prop.bs(.05, f, n)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(3, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.agree2 returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL",  "UL"
  )
  
  res <- ci.agree2(.05, 75, 70, 60, 45, 2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(3, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("power.prop returns valid matrix", {
  colnames_expected <- c(
    "Power"
  )
  
  res <- power.prop(.05, 40, .5, .2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("power.prop2 returns valid matrix", {
  colnames_expected <- c(
    "Power"
  )
  
  res <- power.prop2(.05, 60, 40, .5, .5, .2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("power.prop.ps returns valid matrix", {
  colnames_expected <- c(
    "Power"
  )
  
  res <- power.prop.ps(.05, 45, .5, .5, .4, .2)

  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.prop.inv returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE",        "LL",        "UL"
  )
  
  res <- ci.prop.inv(.05, 5, 67)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.prop2.inv returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE",        "LL",        "UL"
  )
  
  res <- ci.prop2.inv(.05, 10, 10, 48, 213)

  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.agree.3rater returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "LL",        "UL"
  )
  
  f <- c(100, 6, 4, 40, 20, 1, 9, 120)
  res <- ci.agree.3rater(.05, f)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(7, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})



test_that("ci.bayes.prop returns valid matrix", {
  colnames_expected <- c(
    "Posterior mean", "Posterior SD", "LL",        "UL"
  )
  
  res <- ci.bayes.prop(.05, .4, .1, 12, 100)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.pv returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "LL", "UL"
  )
  
  res <- ci.pv(.05, 89, 5, 100, 100, .16)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(2, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})

test_that("ci.prop.fpc returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL",        "UL"
  )
  
  res <- ci.prop.fpc(.05, 12, 100, 400)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.poisson returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "SE", "LL",        "UL"
  )
  
  res <- ci.poisson(.05, 23, 5.25)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.ratio.poisson2 returns valid matrix", {
  colnames_expected <- c(
    "Estimate", "LL",        "UL"
  )
  
  res <- ci.ratio.poisson2(.05, 19, 5, 30, 40.5)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})




test_that("pi.prop returns valid matrix", {
  colnames_expected <- c(
    "LL", "UL"
  )
  
  res <- pi.prop(.1, .225, 80, 120)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("size.ci.tetra returns valid matrix", {
  colnames_expected <- c(
    "Sample size"
  )
  
  res <- size.ci.tetra(.05, .4, .3, .5, .3)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(res[[1,1]], 296)
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("size.ci.prop.prior returns valid matrix", {
  colnames_expected <- c(
    "Sample size"
  )
  
  res <- size.ci.prop.prior(.05, .20, .1425, 200, .1)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(res[[1,1]], 318)
  testthat::expect_equal(colnames(res), colnames_expected)
})

