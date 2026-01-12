# Bayesian credible interval for a semipartial correlation with a skeptical prior

Computes an approximate Bayesian credible interval for a semipartial
correlation with a skeptical prior. The skeptical prior distribution is
Normal with a mean of 0 and a small standard deviation. A skeptical
prior assumes that the population semipartial correlation is within a
range of small values (-r to r). If the skeptic is 95% confident that
the population correlation is between -r and r, then the prior standard
deviation can be set to r/1.96. A semipartial correlation that is less
than .2 in absolute value is typically considered to be "small", and the
prior standard deviation could then be set to .2/1.96. A semipartial
correlation value that is considered to be small will depend on the
application. This function requires the standard error of the estimated
semipartial correlation which can be obtained from the ci.spcor
function.

For more details, see Section 2.36 of Bonett (2021, Volume 2)

## Usage

``` r
ci.bayes.spcor(alpha, prior_sd, cor, se)
```

## Arguments

- alpha:

  alpha level for 1-alpha credibility interval

- prior_sd:

  standard deviation of skeptical prior distribution

- cor:

  estimated semipartial partial correlation

- se:

  standard error of estimated semipartial correlation

## Value

Returns a 1-row matrix. The columns are:

- Posterior mean - posterior mean (Bayesian estimate of correlation)

- LL - lower limit of the credible interval

- UL - upper limit of the credible interval

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
ci.bayes.spcor(.05, .1, .582, .137)
#>  Posterior mean     LL    UL
#>          0.2273 0.0729 0.371

# Should return:
# Posterior mean     LL    UL
#         0.2273 0.0729 0.371

```
