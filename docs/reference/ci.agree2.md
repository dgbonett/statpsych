# Confidence interval for G-index difference in a 2-group design

Computes adjusted Wald confidence intervals for the G-index of agreement
within each group and the difference of G-indices.

For more details, see Section 3.5 of Bonett (2021, Volume 3)

## Usage

``` r
ci.agree2(alpha, n1, f1, n2, f2, k)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- n1:

  sample size (objects) in group 1

- f1:

  number of objects rated in agreement in group 1

- n2:

  sample size (objects) in group 2

- f2:

  number of objects rated in agreement in group 2

- k:

  number of rating categories

## Value

Returns a 3-row matrix. The rows are:

- Row 1: G-index for group 1

- Row 2: G-index for group 2

- Row 3: G-index difference

The columns are:

- Estimate - maximum likelihood estimate of G-index and difference

- SE - standard error

- LL - lower limit of adjusted Wald confidence interval

- UL - upper limit of adjusted Wald confidence interval

## References

Bonett DG (2022). “Statistical inference for G-indices of agreement.”
*Journal of Educational and Behavioral Statistics*, **47**(4), 438–458.
[doi:10.3102/10769986221088561](https://doi.org/10.3102/10769986221088561)
.

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
ci.agree2(.05, 75, 70, 60, 45, 2)
#>         Estimate      SE     LL     UL
#> G1        0.8667 0.02880 0.6975 0.9481
#> G2        0.5000 0.05590 0.2523 0.6852
#> G1 - G2   0.3667 0.06289 0.1117 0.6089

# Should return:
#         Estimate      SE     LL     UL
# G1        0.8667 0.02880 0.6975 0.9481
# G2        0.5000 0.05590 0.2523 0.6852
# G1 - G2   0.3667 0.06289 0.1117 0.6089               

```
