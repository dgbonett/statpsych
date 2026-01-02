# Bayesian credible interval for a Pearson or partial correlation with a skeptical prior

Computes an approximate Bayesian credible interval for a Pearson or
partial correlation with a skeptical prior. The skeptical prior
distribution is Normal with a mean of 0 and a small standard deviation.
A skeptical prior assumes that the population correlation is within a
range of small values (-r to r). If the skeptic is 95% confident that
the population correlation is between -r and r, then the prior standard
deviation can be set to r/1.96. A correlation that is less than .2 in
absolute value is typically considered to be "small", and the prior
standard deviation could then be set to .2/1.96. A correlation value
that is considered to be small will depend on the application. Set s = 0
for a Pearson correlation.

For more details, see Section 1.33 of Bonett (2021, Volume 2)

## Usage

``` r
ci.bayes.cor(alpha, prior_sd, cor, s, n)
```

## Arguments

- alpha:

  alpha level for 1-alpha credibility interval

- prior_sd:

  standard deviation of skeptical prior distribution

- cor:

  estimated Pearson or partial correlation

- s:

  number of control variables

- n:

  sample size

## Value

Returns a 1-row matrix. The columns are:

- Posterior mean - posterior mean (Bayesian estimate of correlation)

- LL - lower limit of the credible interval

- UL - upper limit of the credible interval

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
ci.bayes.cor(.05, .1, .536, 0, 50)
#>  Posterior mean    LL     UL
#>          0.1874 0.028 0.3375

# Should return:
# Posterior mean     LL     UL
#         0.1874  0.028 0.3375

ci.bayes.cor(.05, .1, .536, 0, 300)
#>  Posterior mean     LL     UL
#>          0.4195 0.3352 0.4971

# Should return:
# Posterior mean     LL     UL
#         0.4195 0.3352 0.4971

```
