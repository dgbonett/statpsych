# Confidence interval for a mean with a finite population correction

Computes a confidence interval for a population mean with a finite
population correction (fpc) using the estimated mean, estimated standard
deviation, sample size, and population size. This function is useful
when the sample size is not a small fraction of the population size.

For more details, see Section 1.33 of Bonett (2021, Volume 1)

## Usage

``` r
ci.mean.fpc(alpha, m, sd, n, N)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- m:

  estimated mean

- sd:

  estimated standard deviation

- n:

  sample size

- N:

  population size

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimated mean

- SE - standard error with fpc

- LL - lower limit of the confidence interval with fpc

- UL - upper limit of the confidence interval with fpc

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
ci.mean.fpc(.05, 24.5, 3.65, 40, 300)
#>  Estimate        SE       LL       UL
#>      24.5 0.5381631 23.41146 25.58854

# Should return:
# Estimate        SE       LL       UL
#     24.5 0.5381631 23.41146 25.58854
 
```
