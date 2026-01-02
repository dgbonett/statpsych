# Sample size for a Cronbach reliability confidence interval using an reliability estimate from a prior study

Computes the sample size required to estimate a Cronbach reliability
with desired confidence interval precision in applications where an
estimated Cronbach reliability from a prior study is available. The
actual confidence interval width in the planned study will depend on the
value of the estimated reliability in the planned study. An estimated
Cronbach reliability from a prior study can be used to compute a lower
prediction limit for the estimated reliability in the planned study,
which is then used as a planning value in the sample size analysis. The
probability that the prediction interval will have a width that is less
than the desired width in the planned study is approximately 1 - alpha2.

This sample size approach assumes that the population Cronbach
reliability that was estimated in the prior study is very similar to the
population Cronbach reliability that will be estimated in the planned
study. If an estimated Cronbach reliability from a prior study is not
available, the researcher must use expert opinion to guess the value of
the Cronbach reliability that will be observed in the planned study. The
[size.ci.cronbach](https://dgbonett.github.io/statpsych/reference/size.ci.cronbach.md)
function uses a reliability planning value that is based on expert
opinion regarding the likely value of the reliability estimate that will
be observed in the planned study.

## Usage

``` r
size.ci.cronbach.prior(alpha1, alpha2, rel0, n0, r, w)
```

## Arguments

- alpha1:

  alpha level for 1-alpha1 confidence in the planned study

- alpha2:

  alpha level for the 1-alpha2 prediction interval

- rel0:

  estimated reliability in prior study

- n0:

  sample size in prior study

- r:

  number of measurements (items, raters, forms)

- w:

  desired confidence interval width

## Value

Returns the required sample size

## Examples

``` r
size.ci.cronbach.prior(.05, .10, .86, 50, 5, .15)
#>  Sample size
#>           74

# Should return:
# Sample size
#          71

```
