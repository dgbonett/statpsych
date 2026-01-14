# Fisher confidence interval

Computes a Fisher confidence interval for any type of correlation (e.g.,
Pearson, Spearman, Kendall-tau, tetrachoric, phi, partial, semipartial,
etc.) or ordinal association such as gamma, Somers' d, or tau-b. The
correlation could also be between two latent factors obtained from a SEM
analysis (the Fisher CI will be more accurate than the large-sample CI
from a SEM analysis). The standard error can be a traditional standard
error, a bootstrap standard error, or a robust standard error from a SEM
analysis.

## Usage

``` r
ci.fisher(alpha, cor, se)
```

## Arguments

- alpha:

  alpha value for 1-alpha confidence

- cor:

  estimated correlation or association coefficient

- se:

  standard error of correlation or association coefficient

## Value

Returns a 1-row matrix. The columns are:

- Estimate - correlation (from input)

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## Examples

``` r
ci.fisher(.05, .692, .049)
#>      Estimate     LL     UL
#> [1,]    0.692 0.5833 0.7763

# Should return:
# Estimate     LL     UL
#    0.692 0.5833 0.7763

```
