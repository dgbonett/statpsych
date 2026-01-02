# Confidence intervals for parameters of one-way random effects ANOVA

Computes estimates and confidence intervals for four parameters of the
one-way random effects ANOVA: 1) the superpopulation grand mean, 2) the
square-root within-group variance component, 3) the square-root
between-group variance component, and 4) the omega-squared coefficient.
This function assumes equal sample sizes.

For more details, see Section 3.18 of Bonett (2021, Volume 1)

## Usage

``` r
ci.random.anova(alpha, m, sd, n)
```

## Arguments

- alpha:

  1 - alpha confidence

- m:

  vector of estimated group means

- sd:

  vector of estimated group standard deviations

- n:

  common sample size in each group

## Value

Returns a 4-row matrix. The rows are:

- Grand mean - the mean of the superpopulation of means

- Within SD - the square-root within-group variance component

- Between SD - the square-root between-group variance component

- Omega-squared - the omega-squared coefficient

The columns are:

- Estimate - estimate of parameter

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
m <- c(56.1, 51.2, 60.3, 68.2, 48.9, 70.5)
sd <- c(9.45, 8.79, 9.71, 8.90, 8.31, 9.75)
ci.random.anova(.05, m, sd, 20)
#>                 Estimate         LL         UL
#> Grand mean     59.200000 49.9363896 68.4636104
#> Within SD:      9.166782  8.0509046 10.4373219
#> Between SD:     8.585948  8.3239359  8.8562078
#> Omega-squared:  0.467317  0.2284142  0.8480383

# Should return:
#                 Estimate         LL         UL
# Grand mean     59.200000 49.9363896 68.4636104
# Within SD:      9.166782  8.0509046 10.4373219
# Between SD:     8.585948  8.3239359  8.8562078
# Omega-squared:  0.467317  0.2284142  0.8480383

```
