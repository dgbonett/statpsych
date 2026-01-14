# Sample size for a slope hypothesis test in a general statistical model

Computes the sample size required to test a null hypothesis with desired
power that a population slope coefficient in any general statistical
model is equal to zero. This function requires a standard error estimate
for the slope of interest from a prior or pilot study and the sample
size that was used in the prior or pilot study. This function can be
used for both unstandardized and standardized slopes. This function also
can be used for both unstandardized and standardized factor loadings in
a confirmatory factor analysis model. This function will soon be
replaced with size.test.gen.

## Usage

``` r
size.test.slope.gen(alpha, pow, se, n0, b)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- pow:

  desired power

- se:

  standard error of slope from prior/pilot study

- n0:

  sample size used in prior/pilot study

- b:

  planning value of population slope

## Value

Returns the required sample size

## Examples

``` r
size.test.slope.gen(.05, .8, 3.15, 50, 5)
#>  Sample size
#>          156

# Should return:
#  Sample size
#          156
 
```
