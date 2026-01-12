# Confidence intervals for positive and negative predictive values with retrospective sampling

Computes adjusted Wald confidence intervals for positive and negative
predictive values (PPV and NPV) of a diagnostic test with retrospective
sampling where the population prevalence rate is assumed to be known.
With retrospective sampling, one random sample is obtained from a
subpopulation that is known to have a "positive" outcome, a second
random sample is obtained from a subpopulation that is known to have a
"negative" outcome, and then the diagnostic test (scored "pass" or
"fail") is given in each sample. PPV and NPV can be expressed as a
function of proportion ratios and the known population prevalence rate
(the population proportion who would "pass"). The confidence intervals
for PPV and NPV are based on the Price-Bonett adjusted Wald confidence
interval for a proportion ratio.

For more details, see Section 3.6 of Bonett (2021, Volume 3)

## Usage

``` r
ci.pv(alpha, f1, f2, n1, n2, prev)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- f1:

  number of participants with a positive outcome who pass the test

- f2:

  number of participants with a negative outcome who fail the test

- n1:

  sample size for the positive outcome group

- n2:

  sample size for the negative outcome group

- prev:

  known population proportion with a positive outcome

## Value

Returns a 2-row matrix. The columns are:

- Estimate - adjusted estimate of the predictive value

- LL - lower limit of the adjusted Wald confidence interval

- UL - upper limit of the adjusted Wald confidence interval

## References

Price RM, Bonett DG (2008). “Confidence intervals for a ratio of two
independent binomial proportions.” *Statistics in Medicine*, **27**(26),
5497–5508. ISSN 02776715,
[doi:10.1002/sim.3376](https://doi.org/10.1002/sim.3376) .

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
ci.pv(.05, 89, 5, 100, 100, .16)
#>        Estimate        LL        UL
#> PPV:  0.7640449 0.5838940 0.8819671
#> NPV:  0.9779978 0.9623406 0.9872318

# Should return:
#        Estimate        LL        UL
# PPV:  0.7640449 0.5838940 0.8819671
# NPV:  0.9779978 0.9623406 0.9872318

```
