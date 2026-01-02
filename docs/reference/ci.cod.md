# Confidence interval for a coefficient of dispersion

Computes a confidence interval for a population coefficient of
dispersion (COD). The COD is a mean absolute deviation from the median
divided by the median. The coefficient of dispersion assumes ratio-scale
scores and is a robust alternative to the coefficient of variation (see
[ci.cv](https://dgbonett.github.io/statpsych/reference/ci.cv.md)). An
approximate standard error is recovered from the confidence interval.

For more details, see Section 1.27 of Bonett (2021, Volume 1)

## Usage

``` r
ci.cod(alpha, y)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- y:

  vector of scores

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimated coefficient of dispersion

- SE - recovered standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG, Seier E (2006). “Confidence interval for a coefficient of
dispersion in nonnormal distributions.” *Biometrical Journal*,
**48**(1), 144–148. ISSN 0323-3847,
[doi:10.1002/bimj.200410148](https://doi.org/10.1002/bimj.200410148) .
Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
y <- c(30, 20, 15, 10, 10, 60, 20, 25, 20, 30, 10, 5, 50, 40,
       20, 10, 0, 20, 50)
ci.cod(.05, y)
#>   Estimate        SE        LL       UL
#>  0.5921053 0.1814708 0.3813259 1.092679

# Should return:
#  Estimate        SE        LL       UL
# 0.5921053 0.1814708 0.3813259 1.092679

```
