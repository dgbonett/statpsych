# Sample size for a 2-group Cronbach reliability difference confidence interval

Computes the sample size per group (assuming equal sample sizes)
required to estimate a difference in population Cronbach reliability
coefficients with desired precision in a 2-group design.

## Usage

``` r
size.ci.cronbach2(alpha, rel1, rel2, r, w)
```

## Arguments

- alpha:

  alpha level for hypothesis test

- rel1:

  reliability planning value for group 1

- rel2:

  reliability planning value for group 2

- r:

  number of measurements (items, raters, forms)

- w:

  desired confidence interval width

## Value

Returns the required sample size for each group

## Details

For more details, see Section 2.19 of Bonett (2021, Volume 4)

## References

Bonett DG, Wright TA (2015). “Cronbach's alpha reliability: Interval
estimation, hypothesis testing, and sample size planning.” *Journal of
Organizational Behavior*, **36**(1), 3–15. ISSN 08943796,
[doi:10.1002/job.1960](https://doi.org/10.1002/job.1960) . Bonett DG
(2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
size.ci.cronbach2(.05, .85, .70, 8, .15)
#>  Sample size per group
#>                    180

# Should return:
# Sample size per group
#                   180
 
```
