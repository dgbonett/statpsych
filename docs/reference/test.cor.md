# Hypothesis test for a Pearson or partial correlation

Computes a t test for a test of the null hypothesis that a population
Pearson or partial correlations is equal to 0, or a z test using a
Fisher transformation for a test of the null hypothesis that a Pearson
or partial correlation is equal to some specified nonzero value. Set s =
0 for a Pearson correlation. The hypothesis testing results should be
accompanied with a confidence interval for the population Pearson or
partial correlation value (see
[ci.cor](https://dgbonett.github.io/statpsych/reference/ci.cor.md)).

For more details, see Section 1.19 of Bonett (2021, Volume 2)

## Usage

``` r
test.cor(cor, n, s, h)
```

## Arguments

- cor:

  estimated correlation

- n:

  sample size

- s:

  number of control variables

- h:

  null hypothesis value of correlation

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimate of correlation

- t or z - t test statistic (for h = 0) or z test statistic (for nonzero
  h)

- p - two-sided p-value

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
test.cor(.484, 100, 0, .2)
#>  Estimate      z p
#>     0.484 3.2054 0

# Should return:
# Estimate      z       p
#    0.484 3.2054 0.00135


test.cor(.372, 100, 0, 0)
#>  Estimate      t df       p
#>     0.372 3.9673 98 0.00014

# Should return:
#  Estimate      t df       p
#     0.372 3.9673 98 0.00014

```
