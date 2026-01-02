# Confidence interval for a proportion with a finite population correction

Computes an adjusted Wald interval for a population proportion with a
finite population correction (fpc). This confidence interval is useful
when the sample size is not a small fraction of the population size.

For more details, see Section 1.20 of Bonett (2021, Volume 3)

## Usage

``` r
ci.prop.fpc(alpha, f, n, N)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- f:

  number of participants who have the attribute

- n:

  sample size

- N:

  population size

## Value

Returns a 1-row matrix. The columns are:

- Estimate - adjusted estimate of proportion

- SE - adjusted standard error with fpc

- LL - lower limit of the confidence interval with fpc

- UL - upper limit of the confidence interval with fpc

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
ci.prop.fpc(.05, 61, 75, 415)
#>   Estimate         SE       LL        UL
#>  0.7974684 0.04097594 0.717157 0.8777797

# Should return:
#   Estimate         SE       LL        UL
#  0.7974684 0.04097594 0.717157 0.8777797

```
