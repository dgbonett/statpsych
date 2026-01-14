# Sample size for an intraclass correlation confidence interval using a planning value from a prior study

Computes the sample size required to estimate an intraclass correlation
with desired confidence interval precision in applications where an
estimated intraclass correlation from a prior study is available. The
actual confidence interval width in the planned study will depend on the
value of the estimated intraclas correlation in the planned study. An
estimated intraclass correlation from a prior study can be used to
compute a lower prediction limit for the estimated intraclass
correlation in the planned study, which is then used as a planning value
in the sample size analysis. The probability that the prediction
interval will have a width that is less than the desired width in the
planned study is approximately 1 - alpha2.

This sample size approach assumes that the population intraclass
correlation that was estimated in the prior study is very similar to the
population intraclass correlation that will be estimated in the planned
study. If an estimated intraclas correlation from a prior study is not
available, the researcher must use expert opinion to guess the value of
the intraclass correlation that will be observed in the planned study.
The
[size.ci.icc](https://dgbonett.github.io/statpsych/reference/size.ci.icc.md)
function uses an intraclass correlation planning value that is based on
expert opinion regarding the likely value of the intraclass correlation
estimate that will be observed in the planned study.

## Usage

``` r
size.ci.icc.prior(alpha1, alpha2, cor0, n0, r, w)
```

## Arguments

- alpha1:

  alpha level for 1-alpha1 confidence in the planned study

- alpha2:

  alpha level for the 1-alpha2 prediction interval

- cor0:

  estimated correlation in prior study

- n0:

  sample size in prior study

- r:

  number of measurements (raters, forms)

- w:

  desired confidence interval width

## Value

Returns the required sample size

## Examples

``` r
size.ci.icc.prior(.05, .10, .674, 50, 3, .2)
#>  Sample size
#>          114

# Should return:
# Sample size
#         114

```
