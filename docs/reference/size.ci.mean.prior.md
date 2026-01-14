# Sample size for a mean confidence interval using an estimated variance from a prior study

Computes the sample size required to estimate a population mean with
desired confidence interval precision in applications where an estimated
variance from a prior study is available. The actual confidence interval
width in the planned study will depend on the value of the estimated
variance in the planned study. An estimated variance from a prior study
can be used to compute an upper prediction limit for the estimated
variance in the planned study. The upper prediction limit is then used
as the variance planning value. The probability that the prediction
interval in the planned study will have a width that is less than the
desired width is approximately 1 - alpha2.

This sample size approach assumes that the population variance in the
prior study is very similar to the population variance in the planned
study. If an estimated variance from a prior study is not available, the
researcher must use expert opinion to guess the value of the variance
that will be observed in the planned study. The
[size.ci.mean](https://dgbonett.github.io/statpsych/reference/size.ci.mean.md)
function uses a variance planning value that is based on expert opinion
regarding the likely value of the variance estimate that will be
observed in the planned study.

For more details, see Section 1.31 of Bonett (2021, Volume 1)

## Usage

``` r
size.ci.mean.prior(alpha1, alpha2, var0, n0, w)
```

## Arguments

- alpha1:

  alpha level for 1-alpha1 confidence in the planned study

- alpha2:

  alpha level for the 1-alpha2 prediction interval

- var0:

  estimated variance in prior study

- n0:

  sample size in prior study

- w:

  desired confidence interval width

## Value

Returns the required sample size

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
size.ci.mean.prior(.05, .10, 0.71, 204, .4)
#>  Sample size
#>           88

# Should return:
# Sample size
#          88
```
