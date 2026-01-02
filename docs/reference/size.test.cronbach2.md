# Sample size to test equality of Cronbach reliability coefficients in a 2-group design

Computes the sample size required to test a difference in population
Cronbach reliability coefficients with desired power in a 2-group
design.

## Usage

``` r
size.test.cronbach2(alpha, pow, rel1, rel2, r)
```

## Arguments

- alpha:

  alpha level for hypothesis test

- pow:

  desired power

- rel1:

  reliability planning value for group 1

- rel2:

  reliability planning value for group 2

- r:

  number of measurements (items, raters, forms)

## Value

Returns the required sample size for each group

## References

Bonett DG, Wright TA (2015). “Cronbach's alpha reliability: Interval
estimation, hypothesis testing, and sample size planning.” *Journal of
Organizational Behavior*, **36**(1), 3–15. ISSN 08943796,
[doi:10.1002/job.1960](https://doi.org/10.1002/job.1960) .

## Examples

``` r
size.test.cronbach2(.05, .80, .85, .70, 8)
#>  Sample size per group
#>                     77

# Should return:
# Sample size per group
#                    77
 
```
