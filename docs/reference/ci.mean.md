# Confidence interval for a mean

Computes a confidence interval for a population mean using the estimated
mean, estimated standard deviation, and sample size. Use the t.test
function for raw data input.

For more details, see Section 1.7 of Bonett (2021, Volume 1)

## Usage

``` r
ci.mean(alpha, m, sd, n)
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

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimated mean

- SE - standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Snedecor GW, Cochran WG (1989). *Statistical Methods*, 8th edition. ISU
University Pres, Ames, Iowa.

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
ci.mean(.05, 38.3, 8.14, 10)
#>  Estimate       SE       LL       UL
#>      38.3 2.574094 32.47699 44.12301

# Should return:
# Estimate        SE       LL       UL
#     38.3  2.574094 32.47699 44.12301
 
```
