# Sample size for a paired-samples mean difference confidence interval using an estimated variance and correlation from a prior study

Computes the sample size required to estimate a population mean
difference in a paired-samples design with desired confidence interval
precision in applications where an estimated variance and correlation
from a prior study is available. The actual confidence interval width in
the planned study will depend on the value of the estimated variance of
the difference scores in the planned study. An estimated variance and
correlation from a prior study can be used to compute an upper
prediction limit for the estimated difference score variance in the
planned study. The upper prediction limit is then used as the difference
score variance planning value. The probability that the 1 - alpha1
confidence interval in the planned study will have a width that is less
than the desired width is approximately 1 - alpha2 where alpha1 and
alpha2 are specified values.

This sample size approach assumes that the population variance and
correlation in the prior study is very similar to the population
variance and correlation in the planned study. If information from a
prior study is not available, the researcher must use expert opinion to
guess the values of the variance and correlation that will be observed
in the planned study. The
[size.ci.mean.ps](https://dgbonett.github.io/statpsych/reference/size.ci.mean.ps.md)
function uses variance and correlation planning values that are based on
expert opinion regarding the likely values of the variance and
correlation estimates that will be observed in the planned study.

For more details, see Section 1.31 of Bonett (2021, Volume 1)

## Usage

``` r
size.ci.mean.ps.prior(alpha1, alpha2, var0, cor0, n0, w)
```

## Arguments

- alpha1:

  alpha level for 1-alpha1 confidence in the planned study

- alpha2:

  alpha level for the 1-alpha2 prediction interval

- var0:

  estimated variance in prior study

- cor0:

  estimated correlation in prior study

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
size.ci.mean.ps.prior(.05, .10, 15.2, .78, 10, 2)
#>  Sample size
#>           59

# Should return:
# Sample size
#          59
```
