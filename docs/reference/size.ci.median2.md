# Sample size for a 2-group median difference confidence interval

Computes the sample size for each group required to estimate a
population median difference with desired confidence interval precision
in a 2-group design. Set the variance planning value to the largest
value within a plausible range for a conservatively large sample size.
The sample size requirement depends on the shape of the distribution.
Select one of the four distribution options (Normal, Logistic, Laplace,
Exponential) that approximates the most likely distribution shape in the
planned study. Select the Normal distribution for a conservatively large
sample size requirement. Set R = 1 for equal sample sizes.

For more details, see Sections 1.28 and 2.13 of Bonett (2021, Volume 1)

## Usage

``` r
size.ci.median2(alpha, var, w, R, dist)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- var:

  planning value of average within-group variance

- w:

  desired confidence interval width

- R:

  n2/n1 ratio

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
. Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
size.ci.median2(.05, 37.1, 5, 1, 1)
#>  n1 n2
#>  72 72

# Should return:
# n1  n2
# 72  72

size.ci.median2(.05, 37.1, 5, 2, 4)
#>  n1  n2
#>  51 102

# Should return:
# n1  n2
# 51 102
```
