# Sample size for a median confidence interval

Computes the sample size required to estimate a population median with
desired confidence interval precision. Set the variance planning value
to the largest value within a plausible range for a conservatively large
sample size. The sample size requirement depends on the shape of the
distribution. Select one of the four distribution options (Normal,
Logistic, Laplace, Exponential) that approximates the most likely
distribution shape in the planned study. Select the Normal distribution
for a conservatively large sample size requirement.

For more details, see Section 1.28 of Bonett (2021, Volume 1)

## Usage

``` r
size.ci.median(alpha, var, w, dist)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- var:

  planning value of response variable variance

- w:

  desired confidence interval width

- dist:

  - set to 1 for Normal distribution (skew = 0, kurtosis = 3)

  - set to 2 for Logistic distribution (skew = 0, kurtosis = 4.2)

  - set to 3 for Laplace distribution (skew = 0, kurtosis = 6)

  - set to 4 for Gamma(5) (skew = .89, kurtosis = 4.2)

  - set to 5 for Exponential distribution (skew = 2, kurtosis = 9)

## Value

Returns the required sample size

## References

Bonett DG, Price RM (2002). “Statistical inference for a linear function
of medians: Confidence intervals, hypothesis testing, and sample size
requirements.” *Psychological Methods*, **7**(3), 370–383. ISSN
1939-1463,
[doi:10.1037/1082-989X.7.3.370](https://doi.org/10.1037/1082-989X.7.3.370)
. Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
size.ci.median(.05, 264.4, 10, 1)
#>  Sample size
#>           64

# Should return:
# Sample size
#          64
 
size.ci.median(.05, 264.4, 10, 3)
#>  Sample size
#>           21

# Should return:
# Sample size
#          21

```
