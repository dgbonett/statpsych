# Confidence intervals for a proportion

Computes adjusted Wald (Agresi-Coull), Wilson, and exact confidence
intervals for a population proportion. The Wilson confidence interval
uses a continuity correction.

For more details, see Section 1.5 of Bonett (2021, Volume 3)

## Usage

``` r
ci.prop(alpha, f, n)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- f:

  number of participants who have the attribute

- n:

  sample size

## Value

Returns a 2-row matrix. The columns of row 1 are:

- Estimate - adjusted estimate of proportion

- SE - standard error of adjusted estimate

- LL - lower limit of the adjusted Wald confidence interval

- UL - upper limit of the adjusted Wald confidence interval

The columns of row 2 are:

- Estimate - ML estimate of proportion

- SE - standard error of ML estimate

- LL - lower limit of the Wilson confidence interval

- UL - upper limit of the Wilson confidence interval

The columns of row 3 are:

- Estimate - ML estimate of proportion

- SE - standard error of ML estimate

- LL - lower limit of the exact confidence interval

- UL - upper limit of the exact confidence interval

## References

Agresti A, Coull BA (1998). “Approximate is better than 'exact' for
interval estimation of binomial proportions.” *The American
Statistician*, **52**(2), 119–126. ISSN 0003-1305,
[doi:10.1080/00031305.1998.10480550](https://doi.org/10.1080/00031305.1998.10480550)
. Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
ci.prop(.05, 120, 300)
#>                 Estimate         SE        LL        UL
#> Adjusted Wald  0.4013158 0.02811287 0.3462156 0.4564160
#> Wilson with cc 0.4000000 0.02828427 0.3445577 0.4580464
#> Exact          0.4000000 0.02828427 0.3441290 0.4578664

# Should return:
#                 Estimate         SE        LL        UL
# Adjusted Wald  0.4013158 0.02811287 0.3462156 0.4564160
# Wilson with cc 0.4000000 0.02828427 0.3445577 0.4580464
# Exact          0.4000000 0.02828427 0.3441290 0.4578664

```
