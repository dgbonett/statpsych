# Bayesian credible interval for a normal prior distribution

Computes an approximate Bayesian credible interval for a normal prior
distribution. This function can be used with any parameter estimator
(e.g., mean, mean difference, linear contrast of means, slope
coefficient, standardized mean difference, standardized linear contrast
of means, median, median difference, linear contrast of medians, etc.)
that has an approximate normal sampling distribution. The mean and
standard deviation of the posterior normal distribution are also
reported.

For more details, see Section 1.32 of Bonett (2021, Volume 1)

## Usage

``` r
ci.bayes.normal(alpha, prior_mean, prior_sd, est, se)
```

## Arguments

- alpha:

  alpha level for 1-alpha credibility interval

- prior_mean:

  mean of prior Normal distribution

- prior_sd:

  standard deviation of prior Normal distribution

- est:

  sample estimate

- se:

  standard error of sample estimate

## Value

Returns a 1-row matrix. The columns are:

- Posterior mean - posterior mean of Normal distribution

- Posterior SD - posterior standard deviation of Normal distribution

- LL - lower limit of the credible interval

- UL - upper limit of the credible interval

## References

Gelman A, Carlin JB, Stern HS, Rubin DB (2004). *Bayesian Data
Analysis*, 2nd edition. Chapman & Hall.

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
ci.bayes.normal(.05, 50, 5, 38.3, 2.57)
#>  Posterior mean Posterior SD       LL       UL
#>        40.74511     2.285735 36.26515 45.22506

# Should return:
# Posterior mean Posterior SD       LL       UL
#       40.74511     2.285735 36.26515 45.22506

```
