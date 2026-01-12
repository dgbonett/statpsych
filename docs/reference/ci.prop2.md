# Confidence interval for a 2-group proportion difference

Computes an adjusted Wald confidence interval for a population
proportion difference in a 2-group design.

For more details, see Section 2.2 of Bonett (2021, Volume 3)

## Usage

``` r
ci.prop2(alpha, f1, f2, n1, n2)
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

- Estimate - adjusted estimate of proportion difference

- SE - adjusted standard error

- LL - lower limit of the adjusted Wald confidence interval

- UL - upper limit of the adjusted Wald confidence interval

## References

Agresti A, Caffo B (2000). “Simple and effective confidence intervals
for proportions and differences of proportions result from adding two
successes and two failures.” *The American Statistician*, **54**(4),
280-288. ISSN 00031305,
[doi:10.2307/2685779](https://doi.org/10.2307/2685779) .

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
ci.prop2(.05, 57, 15, 100, 100)
#>   Estimate         SE        LL        UL
#>  0.4117647 0.06083948 0.2925215 0.5310079

# Should return:
#   Estimate         SE        LL        UL
#  0.4117647 0.06083948 0.2925215 0.5310079

```
