# Confidence interval for a biserial-phi correlation

Computes a confidence interval for a population biserial-phi correlation
using a transformation of a confidence interval for an odds ratio with
.5 added to each cell frequency. This measure of association assumes the
group variable is naturally dichotomous and the response variable is
artificially dichotomous.

For more details, see Section 3.4 of Bonett (2021, Volume 3)

## Usage

``` r
ci.biphi(alpha, f1, f2, n1, n2)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- f1:

  number of participants in group 1 who have the attribute

- f2:

  number of participants in group 2 who have the attribute

- n1:

  sample size for group 1

- n2:

  sample size for group 2

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimate of biserial-phi correlation

- SE - standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Ulrich R, Wirtz M (2004). “On the correlation of a naturally and an
artificially dichotomized variable.” *British Journal of Mathematical
and Statistical Psychology*, **57**(2), 235–251. ISSN 00071102,
[doi:10.1348/0007110042307203](https://doi.org/10.1348/0007110042307203)
.

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
ci.biphi(.05, 34, 22, 50, 50)
#>  Estimate     SE    LL    UL
#>     0.275 0.1075 0.049 0.464

# Should return:
#  Estimate     SE    LL    UL
#     0.275 0.1075 0.049 0.464

```
