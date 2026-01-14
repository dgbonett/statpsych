# Sample size for a between-subjects mean linear contrast confidence interval

Computes the sample size in each group (assuming equal sample sizes)
required to estimate a linear contrast of population means with desired
confidence interval precision in a between-subjects design. Set the
variance planning value to the largest value within a plausible range
for a conservatively large sample size.

For more details, see Section 3.24 of Bonett (2021, Volume 1)

## Usage

``` r
size.ci.lc.mean.bs(alpha, var, w, v)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- var:

  planning value of average within-group variance

- w:

  desired confidence interval width

- v:

  vector of between-subjects contrast coefficients

## Value

Returns the required sample size for each group

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
v <- c(.5, .5, -.5, -.5)
size.ci.lc.mean.bs(.05, 8.0, 3.0, v)
#>  Sample size per group
#>                     15

# Should return:
# Sample size per group
#                    15
 
```
