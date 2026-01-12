# Sample size for a within-subjects mean linear contrast confidence interval

Computes the sample size required to estimate a linear contrast of
population means with desired confidence interval precision in a
within-subjects design. Set the variance planning value to the largest
value within a plausible range for a conservatively large sample size.
Set the Pearson correlation planning value to the smallest value within
a plausible range for a conservatively large sample size.

For more details, see Section 4.26 of Bonett (2021, Volume 1)

## Usage

``` r
size.ci.lc.mean.ws(alpha, var, cor, w, q)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- var:

  planning value of average variance of the measurements

- cor:

  planning value of average correlation between measurements

- w:

  desired confidence interval width

- q:

  vector of within-subjects contrast coefficients

## Value

Returns the required sample size

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
q <- c(.5, .5, -.5, -.5)
size.ci.lc.mean.ws(.05, 161.9, .77, 4, q)
#>  Sample size
#>           38

# Should return:
# Sample size
#          38
 
```
