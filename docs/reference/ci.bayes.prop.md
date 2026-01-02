# Bayesian credible interval for a proportion

Computes a Bayesian credible interval for a population proportion using
the mean and standard deviation of a prior Beta distribution along with
sample information. The mean and standard deviation of the posterior
Beta distribution are also reported. For a noninformative prior, set the
prior mean to .5 and the prior standard deviation to 1/sqrt(12) (which
corresponds to a Beta(1,1) distribution). The prior variance must be
less than m(1 - m) where m is the prior mean.

For more details, see Section 1.18 of Bonett (2021, Volume 3)

## Usage

``` r
ci.bayes.prop(alpha, prior_mean, prior_sd, f, n)
```

## Arguments

- alpha:

  alpha level for 1-alpha credibility interval

- prior_mean:

  mean of prior Beta distribution

- prior_sd:

  standard deviation of prior Beta distribution

- f:

  number of participants who have the attribute

- n:

  sample size

## Value

Returns a 1-row matrix. The columns are:

- Posterior mean - posterior mean of Beta distribution

- Posterior SD - posterior standard deviation of Beta distribution

- LL - lower limit of the credible interval

- UL - upper limit of the credible interval

## References

Gelman A, Carlin JB, Stern HS, Rubin DB (2004). *Bayesian Data
Analysis*, 2nd edition. Chapman & Hall. Bonett DG (2021url).
*Statistical Methods for Psychologists*.

## Examples

``` r
ci.bayes.prop(.05, .25, .1, 120, 300)
#>  Posterior mean Posterior SD        LL        UL
#>       0.3916208   0.02742595 0.3387206 0.4458133

# Should return:
# Posterior mean Posterior SD        LL        UL
#      0.3916208   0.02742595 0.3387206 0.4458133

```
