# Sample size for a mean linear contrast test in an ANCOVA

Computes the sample size for each group (assuming equal sample sizes)
required to test a linear contrast of population means in an ANCOVA
model with desired power. In a nonexperimental design, the sample size
is affected by the magnitude of covariate mean differences across
groups. The covariate mean differences can be approximated by specifying
the largest standardized covariate mean difference across all pairwise
comparisons and for all covariates. In an experiment, this standardized
mean difference is set to 0. Set the error variance planning value to
the largest value within a plausible range for a conservatively large
sample size.

For more details, see Section 2.29 of Bonett (2021, Volume 2)

## Usage

``` r
size.test.lc.ancova(alpha, pow, evar, es, s, d, v)
```

## Arguments

- alpha:

  alpha level for hypothesis test

- pow:

  desired power

- evar:

  planning value of within-group (error) variance

- es:

  planning value of linear contrast

- s:

  number of covariates

- d:

  largest standardized mean difference for all covariates

- v:

  vector of between-subjects contrast coefficients

## Value

Returns the required sample size for each group

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
v <- c(.25, .25, .25, .25, -1)
size.test.lc.ancova(.05, .9, 17.5, 4.0, 2, 0, v)
#>  Sample size per group
#>                     17

# Should return:
# Sample size per group
#                    17
 
```
