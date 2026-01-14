# Confidence interval for a standardized linear contrast of means in a within-subjects design

Computes confidence intervals for two types of population standardized
linear contrast of means (unweighted standardizer and level 1
standardizer) in a within-subjects design. Equality of variances is not
assumed, but the correlations among the repeated measures are assumed to
be approximately equal.

For more details, see Section 4.7 of Bonett (2021, Volume 1)

## Usage

``` r
ci.lc.stdmean.ws(alpha, m, sd, cor, n, q)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- m:

  vector of estimated means for levels of within-subjects factor

- sd:

  vector of estimated standard deviations for levels of within-subjects
  factor

- cor:

  average estimated correlation of all measurement pairs

- n:

  sample size

- q:

  vector of within-subjects contrast coefficients

## Value

Returns a 2-row matrix. The columns are:

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
.

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
m <- c(33.5, 37.9, 38.0, 44.1)
sd <- c(3.84, 3.84, 3.65, 4.98)
q <- c(.5, .5, -.5, -.5)
ci.lc.stdmean.ws(.05, m, sd, .672, 20, q)
#>                          Estimate adj Estimate      SE      LL      UL
#> Unweighted standardizer:  -1.3013      -1.2666 0.31479 -1.9182 -0.6843
#> Level 1 standardizer:     -1.3932      -1.3375 0.36618 -2.1109 -0.6755

# Should return:
#                          Estimate adj Estimate      SE      LL      UL
# Unweighted standardizer:  -1.3013      -1.2666 0.31479 -1.9182 -0.6843
# Level 1 standardizer:     -1.3932      -1.3375 0.36618 -2.1109 -0.6755

```
