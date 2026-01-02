# Confidence interval for a standardized linear contrast of means in a between-subjects design

Computes confidence intervals for a population standardized linear
contrast of means in a between-subjects design. The unweighted
standardizer is recommended in experimental designs. The weighted
standardizer is recommended in nonexperimental designs with simple
random sampling. The group 1 standardizer is useful in both experimental
and nonexperimental designs. Equality of variances is not assumed.

For more details, see Section 3.4 of Bonett (2021, Volume 1)

## Usage

``` r
ci.lc.stdmean.bs(alpha, m, sd, n, v)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- m:

  vector of estimated group means

- sd:

  vector of estimated group standard deviation

- n:

  vector of sample sizes

- v:

  vector of between-subjects contrast coefficients

## Value

Returns a 3-row matrix. The columns are:

- Estimate - estimated standardized linear contrast

- adj Estimate - bias adjusted standardized linear contrast estimate

- SE - standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG (2008). “Confidence intervals for standardized linear
contrasts of means.” *Psychological Methods*, **13**(2), 99–109. ISSN
1939-1463,
[doi:10.1037/1082-989X.13.2.99](https://doi.org/10.1037/1082-989X.13.2.99)
. Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
m <- c(6.94, 7.15, 4.60, 3.68)
sd <- c(2.21, 2.83, 2.29, 1.90)
n <- c(40, 40, 40, 40)
v <- c(.5, .5, -.5, -.5)
ci.lc.stdmean.bs(.05, m, sd, n, v)
#>                          Estimate adj Estimate      SE     LL     UL
#> Unweighted standardizer:   1.2459       1.2399 0.17621 0.9005 1.5912
#> Weighted standardizer:     1.2459       1.2399 0.17313 0.9065 1.5852
#> Group 1 standardizer:      1.3145       1.2890 0.22515 0.8732 1.7558

# Should return:
#                          Estimate adj Estimate      SE     LL     UL
# Unweighted standardizer:   1.2459       1.2399 0.17621 0.9005 1.5912 
# Weighted standardizer:     1.2459       1.2399 0.17313 0.9065 1.5852
# Group 1 standardizer:      1.3145       1.2890 0.22515 0.8732 1.7558

```
