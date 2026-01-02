# Sample size for a linear contrast confidence interval in an ANCOVA

Computes the sample size for each group (assuming equal sample sizes)
required to estimate a population linear contrast of means in an ANCOVA
model with desired confidence interval precision. In a nonexperimental
design, the sample size is affected by the magnitude of covariate mean
differences across groups. The covariate mean differences can be
approximated by specifying the largest standardized covariate mean
difference across all pairwise group differences and for all covariates.
In an experiment, this standardized mean difference should be set to 0.
Set the error variance planning value to the largest value within a
plausible range for a conservatively large sample size.

For more details, see Section 2.28 of Bonett (2021, Volume 2)

## Usage

``` r
size.ci.lc.ancova(alpha, evar, s, d, w, v)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- evar:

  planning value of within group (error) variance

- s:

  number of covariates

- d:

  largest standardized mean difference for all covariates

- w:

  desired confidence interval width

- v:

  vector of between-subjects contrast coefficients

## Value

Returns the required sample size for each group

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
v <- c(.5, .5, -.5, -.5)
size.ci.lc.ancova(.05, 6.4, 1, 0, 3.0, v)
#>  Sample size per group
#>                     13

# Should return:
# Sample size per group
#                    13
 
```
