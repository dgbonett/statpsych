# Sample size for a G-index confidence interval

Computes the sample size required to estimate a population G-index of
agreement for two dichotomous ratings with desired confidence interval
precision. Set the G-index planning value to the smallest value within a
plausible range for a conservatively large sample size.

For more details, see Section 3.12 of Bonett (2021, Volume 3)

## Usage

``` r
size.ci.agree(alpha, G, w)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- G:

  planning value of G-index

- w:

  desired confidence interval width

## Value

Returns the required sample size

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
size.ci.agree(.05, .8, .2)
#>  Sample size
#>          139

# Should return:
# Sample size
#         139

```
