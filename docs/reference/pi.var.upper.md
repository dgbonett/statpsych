# Upper prediction limit for an estimated variance

Computes an upper prediction limit for the estimated variance in a
future study for a planned sample size. The prediction limit uses a
variance estimate from a prior study. Several confidence interval sample
size functions in this package require a planning value of the estimated
variance that is expected in the planned study. The upper variance
prediction limit is useful as a variance planning value for the sample
size required to obtain a confidence interval with desired width. This
strategy for specifying a variance planning value is useful in
applications where the population variance in the prior study is assumed
to be very similar to the population variance in the planned study. This
function will be replaced with pi.var which computes both one-sided and
two-sided prediction limits.

## Usage

``` r
pi.var.upper(alpha, var, n0, n)
```

## Arguments

- alpha:

  alpha value for upper 1-alpha confidence

- var:

  estimated variance from prior study

- n0:

  sample size used to estimate variance

- n:

  planned sample size of future study

## Value

Returns an upper prediction estimate (UL) of an estimated variance in a
future study

## References

Hahn GJ (1972). “Simultaneous prediction intervals to contain the
standard deviations or ranges of future samples from a normal
population.” *Journal of the American Statistical Association*,
**67**(340), 938–942.
[doi:10.1080/01621459.1972.10481322](https://doi.org/10.1080/01621459.1972.10481322)
.

## Examples

``` r
pi.var.upper(.05, 15, 40, 100)
#>       UL
#>  23.9724

# Should return:
#      UL
# 23.9724
 
```
