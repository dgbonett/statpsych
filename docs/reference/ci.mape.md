# Confidence interval for a mean absolute prediction error

Computes a confidence interval for a population mean absolute prediction
error (MAPE) in a general linear model. The MAPE is a more robust
alternative to the residual standard deviation. This function requires a
vector of estimated residuals from a general linear model. This
confidence interval does not assume zero excess kurtosis but does assume
symmetry of the population prediction errors.

For more details, see Section 1.16 of Bonett (2021, Volume 2)

## Usage

``` r
ci.mape(alpha, res, s)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- res:

  vector of residuals

- s:

  number of predictor variables in model

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimated mean absolute prediction error

- SE - standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
res <- c(-2.70, -2.69, -1.32, 1.02, 1.23, -1.46, 2.21, -2.10, 2.56,
      -3.02, -1.55, 1.46, 4.02, 2.34)
ci.mape(.05, res, 1)
#>  Estimate        SE       LL       UL
#>    2.3744 0.3314752 1.806022 3.121654

# Should return:
# Estimate        SE       LL       UL
#   2.3744 0.3314752 1.751678 3.218499
 
```
