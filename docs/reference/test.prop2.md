# Hypothesis test for a 2-group proportion difference

Computes a continuity-corrected z-test for a difference of population
proportions in a 2-group design. A confidence interval for a difference
in population proportions is a recommended supplement to the z-test (see
[ci.prop2](https://dgbonett.github.io/statpsych/reference/ci.prop2.md)).

For more details, see Section 2.4 of Bonett (2021, Volume 3)

## Usage

``` r
test.prop2(f1, f2, n1, n2)
```

## Arguments

- f1:

  number of group 1 participants who have the attribute

- f2:

  number of group 2 participants who have the attribute

- n1:

  sample size for group 1

- n2:

  sample size for group 2

## Value

Returns a 1-row matrix. The columns are:

- Estimate - ML estimate of proportion difference

- z - z test statistic

- p - two-sided p-value

## References

Snedecor GW, Cochran WG (1989). *Statistical Methods*, 8th edition. ISU
University Pres, Ames, Iowa. Bonett DG (2021url). *Statistical Methods
for Psychologists*.

## Examples

``` r
test.prop2(39, 24, 50, 50)
#>  Estimate      z       p
#>       0.3 2.8997 0.00373

# Should return:
# Estimate      z       p
#      0.3 2.8997 0.00373

```
