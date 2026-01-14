# Sample size for a between-subjects median linear contrast confidence interval

Computes the sample size in each group (assuming equal sample sizes)
required to estimate a linear contrast of population medians with
desired confidence interval precision in a between-subjects design. Set
the variance planning value to the largest value within a plausible
range for a conservatively large sample size. The sample size
requirement depends on the shape of the distribution. Select one of the
four distribution options (Normal, Logistic, Laplace, Exponential) that
approximates the most likely distribution shape in the planned study.
Select the Normal distribution for a conservatively large sample size
requirement.

For more details, see Sections 1.28 and 3.24 of Bonett (2021, Volume 1)

## Usage

``` r
size.ci.lc.median.bs(alpha, var, w, v, dist)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- var:

  planning value of average within-group variance

- w:

  desired confidence interval width

- v:

  vector of between-subjects contrast coefficients

- dist:

  - set to 1 for Normal distribution (skew = 0, kurtosis = 3)

  - set to 2 for Logistic distribution (skew = 0, kurtosis = 4.2)

  - set to 3 for Laplace distribution (skew = 0, kurtosis = 6)

  - set to 4 for Gamma(5) (skew = .89, kurtosis = 4.2)

  - set to 5 for Exponential distribution (skew = 2, kurtosis = 9)

## Value

Returns the required sample size for each group

## References

Bonett DG, Price RM (2002). “Statistical inference for a linear function
of medians: Confidence intervals, hypothesis testing, and sample size
requirements.” *Psychological Methods*, **7**(3), 370–383. ISSN
1939-1463,
[doi:10.1037/1082-989X.7.3.370](https://doi.org/10.1037/1082-989X.7.3.370)
.

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
v <- c(.5, .5, -1)
size.ci.lc.median.bs(.05, 5.62, 2.0, v, 1)
#>  Sample size per group
#>                     51

# Should return:
# Sample size per group
#                    51

size.ci.lc.median.bs(.05, 5.62, 2.0, v, 4)
#>  Sample size per group
#>                     48

# Should return:
# Sample size per group
#                    33
 
```
