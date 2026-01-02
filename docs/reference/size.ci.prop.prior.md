# Sample size for a proportion confidence interval using an estimated proportion from a prior study

Computes the sample size required to estimate a population proportion
with desired confidence interval precision in applications where an
estimated proportion from a prior study is available. The actual
confidence interval width in the planned study will depend on the value
of the estimated proportion in the planned study. An estimated
proportion from a prior study is used to predict the value of the
estimated proportion in the planned study, and the predicted proportion
estimate is then used as a planning value in the sample size
computation. The probability that the prediction interval in the planned
study will have a width that is less than the desired width is
approximately 1 - alpha2.

This sample size approach assumes that the population proportion in the
prior study is very similar to the population proportion in the planned
study. If an estimated proportion from a prior study is not available
the researcher must use expert opinion to guess the value of the
proportion that will be observed in the planned study. The
[size.ci.prop](https://dgbonett.github.io/statpsych/reference/size.ci.prop.md))
function uses a proportion planning value that is based on expert
opinion regarding the likely value of the proportion estimate that will
be observed in the planned study.

For more details, see Section 1.16 of Bonett (2021, Volume 3)

## Usage

``` r
size.ci.prop.prior(alpha1, alpha2, p0, n0, w)
```

## Arguments

- alpha1:

  alpha level for 1-alpha1 confidence in the planned study

- alpha2:

  alpha level for the 1-alpha2 prediction interval

- p0:

  estimated proportion in prior study

- n0:

  sample size in prior study

- w:

  desired confidence interval width

## Value

Returns the required sample size

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
size.ci.prop.prior(.05, .10, .78, 50, .1)
#>  Sample size
#>          384

# Should return:
# Sample size
#         384

```
