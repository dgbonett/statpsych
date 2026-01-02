# Sample size for a Pearson correlation confidence interval using an estimated correlation from a prior study

Computes the sample size required to estimate a Pearson correlation with
desired confidence interval precision in applications where an estimated
Pearson correlation from a prior study is available. The actual
confidence interval width in the planned study will depend on the value
of the estimated correlation in the planned study. An estimated
correlation from a prior study can be used to compute a prediction
interval for the value of the estimated correlation in the planned
study, which is then used as a planning value in the sample size
analysis. If the prediction interval includes 0, then the correlation
planning value is set to 0; otherwise, the correlation planning value is
set to the lower prediction limit (if the prior correlation is positive)
or the upper prediction limit (if the prior correlation is negative).
The probability that the prediction interval will have a width that is
less than the desired width in the planned study is approximately 1 -
alpha2.

This sample size approach assumes that the population Pearson
correlation that was estimated in the prior study is very similar to the
population Pearson correlation that will be estimated in the planned
study. If an estimated Pearson correlation from a prior study is not
available the researcher must use expert opinion to guess the value of
the Pearson correlation that will be observed in the planned study. The
[size.ci.cor](https://dgbonett.github.io/statpsych/reference/size.ci.cor.md)
function uses a correlation planning value that is based on expert
opinion regarding the likely value of the correlation estimate that will
be observed in the planned study.

For more details, see Section 1.26 of Bonett (2021, Volume 2)

## Usage

``` r
size.ci.cor.prior(alpha1, alpha2, cor0, n0, w)
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

- w:

  desired confidence interval width

## Value

Returns the required sample size

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
size.ci.cor.prior(.05, .10, -.56, 120, .2)
#>  Sample size
#>          246

# Should return:
# Sample size
#         246

```
