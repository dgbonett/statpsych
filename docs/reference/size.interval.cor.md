# Sample size for an interval test of a Pearson or partial correlation

Computes the sample size required to perform an interval test for a
population Pearson or a partial correlation with desired power where the
interval midpoint is equal to zero. This function can be used to plan a
study where the goal is to show that the population correlation is
small. Set s = 0 for a Pearson correlation. The correlation planning
value must be a value within the hypothesized interval.

For more details, see Section 1.25 of Bonett (2021, Volume 2)

## Usage

``` r
size.interval.cor(alpha, pow, cor, s, h)
```

## Arguments

- alpha:

  alpha level for hypothesis test

- pow:

  desired power

- cor:

  planning value of correlation

- s:

  number of control variables

- h:

  upper limit of hypothesized interval

## Value

Returns the required sample size

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
size.interval.cor(.05, .9, .1, 0, .3)
#>  Sample size
#>          251

# Should return:
# Sample size
#         251
 
```
