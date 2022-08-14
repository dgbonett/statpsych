library(statpsych)
# 
# 
# test_that("ci.cramer returns valid matrix", {
#   colnames_expected <- c(
#     "Cramer's V",     "SE",     "LL",     "UL"
#   )
#   
#   res <- ci.cramer(.05, 19.21, 2, 3, 200)
#   
#   testthat::expect_equal(class(res), c("matrix", "array"))
#   testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
#   testthat::expect_equal(colnames(res), colnames_expected)
# })
# 
# 
# test_that("sim.ci.mean1 returns valid matrix", {
#   colnames_expected <- c(
#      "Coverage", "Lower Error", "Upper Error", "Ave CI Width"
#   )
# 
#   res <- sim.ci.mean1(.05, 40, 4)
# 
#   testthat::expect_equal(class(res), c("matrix", "array"))
#   testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
#   testthat::expect_equal(colnames(res), colnames_expected)
# })
# 
# 
# test_that("sim.ci.mean2 returns valid matrix", {
#   colnames_expected <- c(
#     "Coverage", "Lower Error", "Upper Error", "Ave CI Width"
#   )
# 
#   res <- sim.ci.mean2(.05, 30, 25, 1.5, 4, 5)
# 
#   testthat::expect_equal(class(res), c("matrix", "array"))
#   testthat::expect_equal(dim(res), c(2, length(colnames_expected)))
#   testthat::expect_equal(colnames(res), colnames_expected)
# })
# 
# 
# test_that("sim.ci.mean.ps returns valid matrix", {
#   colnames_expected <- c(
#     "Coverage", "Lower Error", "Upper Error", "Ave CI Width"
#   )
# 
#   res <- sim.ci.mean.ps(.05, 30, 1.5, .7, 4, 5)
# 
#   testthat::expect_equal(class(res), c("matrix", "array"))
#   testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
#   testthat::expect_equal(colnames(res), colnames_expected)
# })
# 
# 
# test_that("sim.ci.median1 returns valid matrix", {
#   colnames_expected <- c(
#     "Coverage", "Lower Error", "Upper Error", "Ave CI Width"
#   )
# 
#   res <- sim.ci.median1(.05, 20, 5)
# 
#   testthat::expect_equal(class(res), c("matrix", "array"))
#   testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
#   testthat::expect_equal(colnames(res), colnames_expected)
# })
# 
# 
# test_that("sim.ci.cor returns valid matrix", {
#   colnames_expected <- c(
#     "Coverage", "Lower Error", "Upper Error", "Ave CI Width"
#   )
# 
#   res <- sim.ci.cor(.05, 30, .7, 4, 5)
# 
#   testthat::expect_equal(class(res), c("matrix", "array"))
#   testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
#   testthat::expect_equal(colnames(res), colnames_expected)
# })
# 
# 
# test_that("sim.ci.spear returns valid matrix", {
#   colnames_expected <- c(
#     "Coverage", "Lower Error", "Upper Error", "Ave CI Width"
#   )
# 
#   res <- sim.ci.spear(.05, 30, .7, 4, 5)
# 
#   testthat::expect_equal(class(res), c("matrix", "array"))
#   testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
#   testthat::expect_equal(colnames(res), colnames_expected)
# })


test_that("ci.2x2.mean.mixed returns valid matrix", {
  colnames_expected <- c(
    "Estimate",        "SE",         "t",       "df",            "p",         "LL",        "UL"
  )
  
  
  y11 = c(18, 19, 20, 17, 20, 16)
  y12 = c(19, 18, 19, 20, 17, 16)
  y21 = c(19, 16, 16, 14, 16, 18)
  y22 = c(16, 10, 12,  9, 13, 15)
  res <- ci.2x2.mean.mixed(.05, y11, y12, y21, y22)

  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(7, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.2x2.mean.ws returns valid matrix", {
  colnames_expected <- c(
    "Estimate",        "SE",         "t",       "df",            "p",         "LL",        "UL"
  )
  
  
  y11 = c(1,2,3,4,5,7,7)
  y12 = c(1,0,2,4,3,8,7)
  y21 = c(4,5,6,7,8,9,8)
  y22 = c(5,6,8,7,8,9,9)
  res <- ci.2x2.mean.ws(.05, y11, y12, y21, y22)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(7, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.2x2.mean.bs returns valid matrix", {
  colnames_expected <- c(
    "Estimate",        "SE",         "t",       "df",            "p",         "LL",        "UL"
  )
  
  
  y11 = c(14, 15, 11, 7, 16, 12, 15, 16, 10, 9)
  y12 = c(18, 24, 14, 18, 22, 21, 16, 17, 14, 13)
  y21 = c(16, 11, 10, 17, 13, 18, 12, 16, 6, 15)
  y22 = c(18, 17, 11, 9, 9, 13, 18, 15, 14, 11)
  res <- ci.2x2.mean.bs(.05, y11, y12, y21, y22)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(7, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.2x2.prop.bs returns valid matrix", {
  colnames_expected <- c(
    "Estimate",        "SE",         "z",       "p",            "LL",         "UL"
  )
  

  f = c(15, 24, 28, 23)
  n = c(50, 50, 50, 50)
  res <- ci.2x2.prop.bs(.05, f, n)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(7, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.2x2.prop.mixed returns valid matrix", {
  colnames_expected <- c(
    "Estimate",        "SE",         "z",       "p",            "LL",         "UL"
  )
  
  
  group1 = c(23, 42, 24, 11)
  group2 = c(26, 27, 13, 34)
  res <- ci.2x2.prop.mixed (.05, group1, group2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(7, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("ci.2x2.prop.mixed returns valid matrix", {
  colnames_expected <- c(
    "Estimate",        "SE",         "z",       "p",            "LL",         "UL"
  )
  
  
  group1 = c(23, 42, 24, 11)
  group2 = c(26, 27, 13, 34)
  res <- ci.2x2.prop.mixed (.05, group1, group2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(7, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})