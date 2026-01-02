# Hypothesis test for a 2-group Pearson or partial correlation difference

Computes a z test for a difference of population Pearson or partial
correlations in a 2-group design. Set s = 0 for a Pearson correlation.
The hypothesis testing results should be accompanied with a confidence
interval for the difference in population correlation values (see
[ci.cor2](https://dgbonett.github.io/statpsych/reference/ci.cor2.md)).

## Usage

``` r
test.cor2(cor1, cor2, n1, n2, s)
```

## Arguments

- cor1:

  estimated correlation for group 1

- cor2:

  estimated correlation for group 2

- n1:

  sample size for group 1

- n2:

  sample size for group 2

- s:

  number of control variables

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimate of correlation difference

- z - z test statistic

- p - two-sided p-value

## Examples

``` r
test.cor2(.684, .437, 100, 125, 0)
#>  Estimate      z       p
#>     0.247 2.7057 0.00682

# Should return:
# Estimate      z       p
#    0.247 2.7057 0.00682

```
