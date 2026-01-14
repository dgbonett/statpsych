# Confidence interval for a G-index of agreement

Computes an adjusted Wald confidence interval for a G-index of agreement
between two polychotomous ratings. This function requires the number of
objects that were given the same rating by both raters. The G-index
corrects for chance agreement. The G-index is a better measure of
agreement than Cohen's kappa, and the confidence interval for the
G-index used here has better small-sample properties than the confidence
interval for Cohen's kappa.

For more details, see Section 3.5 of Bonett (2021, Volume 3)

## Usage

``` r
ci.agree(alpha, n, f, k)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- n:

  sample size

- f:

  number of objects rated in agreement

- k:

  number of rating categories

## Value

Returns a 1-row matrix. The columns are:

- Estimate - maximum likelihood estimate of G-index

- SE - standard error

- LL - lower limit of the adjusted Wald confidence interval

- UL - upper limit of the adjusted Wald confidence interval

## References

Bonett DG (2022). “Statistical inference for G-indices of agreement.”
*Journal of Educational and Behavioral Statistics*, **47**(4), 438–458.
[doi:10.3102/10769986221088561](https://doi.org/10.3102/10769986221088561)
.

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
ci.agree(.05, 250, 214, 3)
#>  Estimate      SE     LL     UL
#>     0.784 0.03331 0.7098 0.8414

# Should return:
#  Estimate      SE     LL     UL
#     0.784 0.03331 0.7098 0.8414

```
