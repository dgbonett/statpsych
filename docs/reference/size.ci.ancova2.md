# Sample size for a 2-group ANCOVA confidence interval

Computes the sample size for each group required to estimate a mean
difference in a 2-group ANCOVA model with desired confidence interval
precision. In a nonexperimental design, the sample size is affected by
the magnitude of covariate mean differences across groups. The covariate
mean differences can be approximated by specifying the largest
standardized covariate mean difference of all covariates. In an
experiment, this standardized mean difference should be set to 0. Set
the error variance planning value to the largest value within a
plausible range for a conservatively large sample size.

For more details, see Section 2.28 of Bonett (2021, Volume 2)

## Usage

``` r
size.ci.ancova2(alpha, evar, s, d, w, R)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- evar:

  planning value of within group (error) variance

- s:

  number of covariates

- d:

  largest standardized mean difference of all covariates

- w:

  desired confidence interval width

- R:

  ratio of n2/n1

## Value

Returns the required sample size for each group

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
size.ci.ancova2(.05, 1.37, 1, 0, 1.5, 1)
#>  n1 n2
#>  21 21

# Should return:
#  n1 n2
#  21 21

size.ci.ancova2(.05, 1.37, 1, 0, 1.5, 2)
#>  n1 n2
#>  16 32

# Should return:
#  n1 n2
#  16 32

size.ci.ancova2(.05, 1.37, 1, .75, 1.5, 1)
#>  n1 n2
#>  24 24

# Should return:
#  n1 n2
#  24 24
 
```
