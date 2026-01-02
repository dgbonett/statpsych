# Confidence interval for an intraclass reliablity coefficient

Computes a confidence interval for a population intraclass reliability
coefficient using mean squared estimates from a two-way ANOVA. This
function will compute point and interval estimates of the ICC(C, 1) and
ICC(C, r) reliability coefficients where ICC(C, 1) is the reliability of
a single measurements (e.g., a single rater or a single form of a test)
and ICC(C, r) is the reliability of a sum or average of r measurements.
ICC(C, r) is the same as Cronbach's reliability coefficient. The
ci.cronbach function uses a point estimate of Cronbach's reliability as
input. The confidence intervals used in this function assume parallel
measurements which implies a compound symmetric covariance matrix of the
r measurements.

## Usage

``` r
ci.icc(alpha, MSr, MSe, r, n)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- MSr:

  mean square for rows

- MSe:

  error mean square

- r:

  number of measurements (items, raters, forms)

- n:

  sample size

## Value

Returns a 2-row matrix. The first row results are for ICC(C, 1) and the
second row results are for ICC(C, r). The columns are:

- Estimate - estimated reliability

- SE - standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

McGraw KO, Wong SP (1996). “Forming inferences about some intraclass
correlation coefficients.” *Psychological Methods*, **1**(1), 30–46.
[doi:10.1037/1082-989X.1.1.30](https://doi.org/10.1037/1082-989X.1.1.30)
. Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
ci.icc(.05, 48.2, 11.3, 5, 30)
#>           Estimate      SE     LL     UL
#> ICC(C, 1)   0.3951 0.09166 0.2311 0.5853
#> ICC(C, r)   0.7656 0.07005 0.6005 0.8759

# Should return:
#           Estimate      SE     LL     UL
# ICC(C, 1)   0.3951 0.09166 0.2311 0.5853
# ICC(C, r)   0.7656 0.07005 0.6005 0.8759
 
```
