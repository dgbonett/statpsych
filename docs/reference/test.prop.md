# Hypothesis test for a proportion

Computes a continuity-corrected z-test for a population proportion in a
1-group design. A confidence interval for a population proportion is a
recommended supplement to the z-test (see
[ci.prop](https://dgbonett.github.io/statpsych/reference/ci.prop.md)).

For more details, see Section 1.7 of Bonett (2021, Volume 3)

## Usage

``` r
test.prop(f, n, h)
```

## Arguments

- f:

  number of participants who have the attribute

- n:

  sample size

- h:

  null hypothesis value of proportion

## Value

Returns a 1-row matrix. The columns are:

- Estimate - ML estimate of proportion

- z - z test statistic

- p - two-sided p-value

## References

Snedecor GW, Cochran WG (1989). *Statistical Methods*, 8th edition. ISU
University Pres, Ames, Iowa.

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
test.prop(76, 100, .6)
#>  Estimate      z       p
#>   0.00156 3.1639 0.00156

# Should return:
# Estimate      z       p
#     0.76 3.1639 0.00156

```
