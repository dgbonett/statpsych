# Prediction limits for sample value of Cronbach reliability in a future study

Computes approximate one-sided or two-sided prediction limits for the
estimated Cronbach reliability in a future study with a planned sample
size of n. The prediction interval uses a Cronbach reliability estimate
from a prior study.

The size.ci.cronbach and size.ci.cronbach2 functions require a planning
value of the expected sample value of Cronbach's reliability in the
planned study. A one-sided lower prediction limit for the sample value
of Cronbach's reliability in the planned study can be used in the
size.ci.cronbach and size.ci.cronbach2 functions to obtain
conservatively large sample size requirements. This strategy for
specifying a reliability planning value is useful in applications where
the population Cronbach reliability in the prior study is assumed to be
very similar to the population Cronbach reliability in the planned
study.

## Usage

``` r
pi.cronbach(alpha, rel, r, n0, n, type)
```

## Arguments

- alpha:

  alpha value for 1-alpha confidence

- rel:

  estimated Cronbach reliability from prior study

- r:

  number of measurements (e.g., items, raters, etc.)

- n0:

  sample size used to estimate reliability in prior study

- n:

  planned sample size of future study

- type:

  - set to 1 for two-sided prediction interval

  - set to 2 for one-sided upper prediction limit

  - set to 3 for one-sided lower prediction limit

## Value

Returns one-sided or two-sided prediction limit(s) of an estimated
Cronbach reliability in a future study

## Examples

``` r
pi.cronbach(.1, .852, 5, 100, 150, 1)
#>      LL     UL
#>  0.7923 0.8945

# Should return:
#       LL     UL
#   0.7944 0.8956
 
pi.cronbach(.1, .852, 5, 100, 150, 3)
#>      LL
#>  0.8073

# Should return:
#      LL
#  0.8092
 
```
