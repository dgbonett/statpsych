# Sample size for a paired-samples sign test

Computes sample size required for a paired-samples sign test with
desired power. The null hypothesis can be expressed in terms of a
population proportion. In a paired-samples experiment, the proportion is
defined as the proportion of members in the population with scores that
would be larger under treatment 1 than treatment 2. In a paired-samples
nonexperiment, the proportion is the proportion of members in the
population with measurement 1 scores that are larger than their
measurement 2 scores. Under the null hypothesis, the population
proportion is equal to to .5. This function requires a planning value of
the population proportion.

## Usage

``` r
size.test.sign.ps(alpha, pow, p)
```

## Arguments

- alpha:

  alpha level for hypothesis test

- pow:

  desired power

- p:

  planning value of proportion

## Value

Returns the required sample size

## Examples

``` r
size.test.sign.ps(.05, .90, .75)
#>  Sample size
#>           42

# Should return:
# Sample size
#          42
 
```
