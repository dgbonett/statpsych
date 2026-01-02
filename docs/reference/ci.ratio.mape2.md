# Confidence interval for a ratio of mean absolute prediction errors in a 2-group design

Computes a confidence interval for a ratio of population mean absolute
prediction errors from a general linear model in two independent groups.
The number of predictor variables can differ across groups and the two
models can be non-nested. This function requires a vector of estimated
residuals from each group. This function does not assume zero excess
kurtosis but does assume symmetry in the population prediction errors
for the two models.

## Usage

``` r
ci.ratio.mape2(alpha, res1, res2, s1, s2)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- res1:

  vector of residuals from group 1

- res2:

  vector of residuals from group 2

- s1:

  number of predictor variables used in group 1

- s2:

  number of predictor variables used in group 2

## Value

Returns a 1-row matrix. The columns are:

- MAPE1 - bias adjusted mean absolute prediction error for group 1

- MAPE2 - bias adjusted mean absolute prediction error for group 2

- MAPE1/MAPE2 - ratio of bias adjusted mean absolute prediction errors

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## Examples

``` r
res1 <- c(-2.70, -2.69, -1.32, 1.02, 1.23, -1.46, 2.21, -2.10, 2.56, -3.02
        -1.55, 1.46, 4.02, 2.34)
res2 <- c(-0.71, -0.89, 0.72, -0.35, 0.33 -0.92, 2.37, 0.51, 0.68, -0.85,
        -0.15, 0.77, -1.52, 0.89, -0.29, -0.23, -0.94, 0.93, -0.31 -0.04)
ci.ratio.mape2(.05, res1, res2, 1, 1)
#>    MAPE1     MAPE2 MAPE1/MAPE2       LL       UL
#>  2.58087 0.8327273    3.099298 1.917003 5.010761

# Should return:
#   MAPE1     MAPE2 MAPE1/MAPE2       LL       UL
# 2.58087 0.8327273    3.099298 1.917003 5.010761
 
```
