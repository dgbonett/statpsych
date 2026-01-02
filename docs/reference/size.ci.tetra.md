# Sample size for a tetrachoric correlation confidence interval

Computes the sample size required to estimate a tetrachoric correlation
with desired confidence interval precision. Set the tetrachoric planning
value to the smallest absolute value within a plausible range for a
conservatively large sample size.

## Usage

``` r
size.ci.tetra(alpha, p1, p2, cor, w)
```

## Arguments

- alpha:

  alpha level for 1 - alpha confidence

- p1:

  planning value for row 1 marginal proportion

- p2:

  planning value for column 1 marginal proportion

- cor:

  tetrachoric planning value

- w:

  desired confidence interval width

## Value

Returns the required sample size

## References

Bonett DG, Price RM (2005). “Inferential methods for the tetrachoric
correlation coefficient.” *Journal of Educational and Behavioral
Statistics*, **30**(2), 213–225. ISSN 1076-9986,
[doi:10.3102/10769986030002213](https://doi.org/10.3102/10769986030002213)
.

## Examples

``` r
size.ci.tetra(.05, .5, .3, .7, .25)
#>  Sample size
#>          304

# Should return:
#  Sample size
#          304

```
