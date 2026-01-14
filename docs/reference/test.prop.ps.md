# Hypothesis test for a paired-samples proportion difference

Computes a continuity-corrected McNemar test for equality of population
proportions in a paired-samples design. This function requires the
frequency counts from a 2 x 2 contingency table for two paired
dichotomous measurements. A confidence interval for a difference in
population proportions (see
[ci.prop.ps](https://dgbonett.github.io/statpsych/reference/ci.prop.ps.md))
is a recommended supplement to the McNemar test.

For more details, see Section 3.2 of Bonett (2021, Volume 3)

## Usage

``` r
test.prop.ps(f00, f01, f10, f11)
```

## Arguments

- f00:

  number participants with y = 0 and x = 0

- f01:

  number participants with y = 0 and x = 1

- f10:

  number participants with y = 1 and x = 0

- f11:

  number participants with y = 1 and x = 1

## Value

Returns a 1-row matrix. The columns are:

- Estimate - ML estimate of proportion difference

- z - z test statistic

- p - two-sided p-value

## References

Snedecor GW, Cochran WG (1989). *Statistical Methods*, 8th edition. ISU
University Pres, Ames, Iowa.

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
test.prop.ps(12, 4, 26, 6)
#>    Estimate      z       p
#>  -0.4583333 3.8341 0.00013

# Should return:
#   Estimate      z       p
# -0.4583333 3.8341 0.00013

```
