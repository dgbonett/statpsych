# Confidence interval for a coefficient of variation

Computes a confidence interval for a population coefficient of variation
(standard deviation divided by mean). This confidence interval is the
reciprocal of a confidence interval for a standardized mean (see
[ci.stdmean](https://dgbonett.github.io/statpsych/reference/ci.stdmean.md)).
An approximate standard error is recovered from the confidence interval.
The coefficient of variation assumes ratio-scale scores.

For more details, see Section 1.27 of Bonett (2021, Volume 1)

## Usage

``` r
ci.cv(alpha, m, sd, n)
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

- Estimate - estimated coefficient of variation

- SE - recovered standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
ci.cv(.05, 24.5, 3.65, 40)
#>   Estimate         SE        LL        UL
#>  0.1489796 0.01817373 0.1214381 0.1926778

# Should return:
#  Estimate        SE        LL       UL
# 0.1489796 0.01817373 0.1214381 0.1926778
 
```
