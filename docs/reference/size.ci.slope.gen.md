# Sample size for a slope confidence interval in a general statistical model

Computes the sample size required to estimate a slope coefficient with
desired confidence interval precision in any type of statistical model.
This function requires a standard error estimate for the slope of
interest from a prior or pilot study and the sample size that was used
in the prior or pilot study. This function can be used for both
unstandardized and standardized slopes. This function also can be used
for both unstandardized and standardized factor loadings in a
confirmatory factor analysis model. This function will soon be replaced
with size.ci.gen.

For more details, see Section 2.28 of Bonett (2021, Volume 2)

## Usage

``` r
size.ci.slope.gen(alpha, se, n0, w)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- se:

  standard error of slope from prior/pilot study

- n0:

  sample size used in prior/pilot study

- w:

  desired confidence interval width

## Value

Returns the required sample size

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
size.ci.slope.gen(.05, 3.15, 50, 5)
#>  Sample size
#>          305

# Should return:
#  Sample size
#          305
 
```
