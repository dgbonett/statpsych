# Sample size for a 2-group ANCOVA hypothesis test

Computes the sample size for each group required to test a mean
difference in an ANCOVA model with desired power in a 2-group design. In
a nonexperimental design, the sample size is affected by the magnitude
of covariate mean differences across groups. The covariate mean
differences can be approximated by specifying the largest standardized
covariate mean difference across of all covariates. In an experiment,
this standardized mean difference is set to 0. Set the error variance
planning value to the largest value within a plausible range for a
conservatively large sample size.

For more details, see Section 2.29 of Bonett (2021, Volume 2)

## Usage

``` r
size.test.ancova2(alpha, pow, evar, es, s, d, R)
```

## Arguments

- alpha:

  alpha level for hypothesis test

- pow:

  desired power

- evar:

  planning value of within-group (error) variance

- es:

  planning value of mean difference

- s:

  number of covariates

- d:

  largest standardized mean difference of all covariates

- R:

  n2/n1 ratio

## Value

Returns the required sample size for each group

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
size.test.ancova2(.05, .9, 1.37, .7, 1, 0, 1)
#>  n1 n2
#>  61 61

# Should return:
#  n1 n2
#  61 61

size.test.ancova2(.05, .9, 1.37, .7, 1, 0, 2)
#>  n1 n2
#>  47 94

# Should return:
#  n1 n2
#  47 94

size.test.ancova2(.05, .9, 1.37, .7, 1, .5, 1)
#>  n1 n2
#>  65 65

# Should return:
#  n1 n2
#  65 65
 
```
