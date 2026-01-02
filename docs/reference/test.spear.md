# Hypothesis test for a Spearman correlation

Computes a t test for a test of the null hypothesis that a population
Spearman correlation is equal to 0, or a z test using a Fisher
transformation for a test of the null hypothesis that a Spearman
correlation is equal to some specified nonzero value. The hypothesis
testing results should be accompanied with a confidence interval for the
population Spearman correlation value (see
[ci.spear](https://dgbonett.github.io/statpsych/reference/ci.spear.md)).

## Usage

``` r
test.spear(cor, h, n)
```

## Arguments

- cor:

  estimated correlation

- h:

  null hypothesis value of correlation

- n:

  sample size

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimate of correlation

- t or z - t test statistic (for h = 0) or z test statistic (for nonzero
  h)

- p - two-sided p-value

## Examples

``` r
test.spear(.471, .2, 100)
#>  Estimate      z p
#>     0.471 3.0096 0

# Should return:
# Estimate       z       p
#     0.471 3.0096 0.00262


test.spear(.341, 0, 100)
#>  Estimate     t df       p
#>     0.341 3.591 98 0.00052

# Should return:
#  Estimate     t df       p
#     0.341 3.591 98 0.00052

```
