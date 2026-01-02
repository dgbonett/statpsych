# Confidence interval for a coefficient of quartile variation

Computes a distribution-free confidence interval for a population
coefficient of quartile variation which is defined as (Q3 - Q1)/(Q3 +
Q1) where Q1 is the 25th percentile and Q3 is the 75th percentile. The
coefficient of quartile variation assumes ratio-scale scores and is a
robust alternative to the coefficient of variation (see
[ci.cv](https://dgbonett.github.io/statpsych/reference/ci.cv.md)). The
25th and 75th percentiles are computed using the type = 2 method (SAS
default).

For more details, see Section 1.27 of Bonett (2021, Volume 1)

## Usage

``` r
ci.cqv(alpha, y)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- y:

  vector of scores

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimated coefficient of quartile variation

- SE - standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG (2006). “Confidence interval for a coefficient of quartile
variation.” *Computational Statistics and Data Analysis*, **50**(11),
2953–2957.
[doi:10.1016/j.csda.2005.05.007](https://doi.org/10.1016/j.csda.2005.05.007)
. Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
y <- c(30, 20, 15, 10, 10, 60, 20, 25, 20, 30, 10, 5, 50, 40,
       20, 10, 0, 20, 50)
ci.cqv(.05, y)
#>  Estimate        SE        LL        UL
#>       0.5 0.1552485 0.2617885 0.8841821

# Should return:
# Estimate        SE        LL       UL
#      0.5 0.1552485 0.2617885 0.8841821

```
