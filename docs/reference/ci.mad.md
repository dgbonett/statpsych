# Confidence interval for a mean absolute deviation

Computes a confidence interval for a population mean absolute deviation
from the median (MAD). The MAD is a robust alternative to the standard
deviation.

For more details, see Section 1.26 of Bonett (2021, Volume 1)

## Usage

``` r
ci.mad(alpha, y)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- y:

  vector of scores

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimated mean absolute deviation

- SE - standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG, Seier E (2003). “Confidence intervals for mean absolute
deviations.” *The American Statistician*, **57**(4), 233–236. ISSN
0003-1305,
[doi:10.1198/0003130032323](https://doi.org/10.1198/0003130032323) .

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
y <- c(30, 20, 15, 10, 10, 60, 20, 25, 20, 30, 10, 5, 50, 40, 
       20, 10, 0, 20, 50)
ci.mad(.05, y)
#>  Estimate       SE       LL       UL
#>      12.5 2.876103 7.962667 19.62282

# Should return:
# Estimate       SE       LL       UL
#     12.5 2.876103 7.962667 19.62282

```
