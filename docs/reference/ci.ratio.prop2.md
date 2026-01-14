# Confidence interval for a 2-group proportion ratio

Computes an adjusted Wald confidence interval for a population
proportion ratio in a 2-group design.

For more details, see Section 2.3 of Bonett (2021, Volume 3)

## Usage

``` r
ci.ratio.prop2(alpha, f1, f2, n1, n2)
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

- Estimate - adjusted estimate of proportion ratio

- LL - lower limit of the adjusted Wald confidence interval

- UL - upper limit of the adjusted Wald confidence interval

## References

Price RM, Bonett DG (2008). “Confidence intervals for a ratio of two
independent binomial proportions.” *Statistics in Medicine*, **27**(26),
5497–5508. ISSN 02776715,
[doi:10.1002/sim.3376](https://doi.org/10.1002/sim.3376) .

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
ci.ratio.prop2(.05, 62, 35, 100, 100)
#>  Estimate       LL       UL
#>  1.765957 1.297164 2.404172

# Should return:
# Estimate       LL       UL
# 1.765957 1.297164 2.404172

```
