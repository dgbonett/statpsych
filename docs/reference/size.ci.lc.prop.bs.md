# Sample size for a between-subjects proportion linear contrast confidence interval

Computes the sample size in each group (assuming equal sample sizes)
required to estimate a linear contrast of population proportions with
desired confidence interval precision in a between-subjects design. Set
the proportion planning values to .5 for a conservatively large sample
size.

For more details, see Section 2.22 of Bonett (2021, Volume 3)

## Usage

``` r
size.ci.lc.prop.bs(alpha, p, w, v)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- p:

  vector of proportion planning values

- w:

  desired confidence interval width

- v:

  vector of between-subjects contrast coefficients

## Value

Returns the required sample size for each group

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
p <- c(.1, .2, .2, .3)
v <- c(.5, -.5, .5, -.5)
size.ci.lc.prop.bs(.05, p, .2, v)
#>  Sample size per group
#>                     60

# Should return:
# Sample size per group
#                    60

```
