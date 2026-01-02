# Confidence interval for a difference in Cronbach reliabilities in a 2-group design

Computes a confidence interval for a difference in population Cronbach
reliability coefficients in a 2-group design. The number of measurements
(items, raters, forms) used in each group need not be equal.

For more details, see Section 2.15 of Bonett (2021, Volume 4)

## Usage

``` r
ci.cronbach2(alpha, rel1, rel2, r1, r2, n1, n2)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- rel1:

  estimated Cronbach reliability for group 1

- rel2:

  estimated Cronbach reliability for group 2

- r1:

  number of measurements used in group 1

- r2:

  number of measurements used in group 2

- n1:

  sample size for group 1

- n2:

  sample size for group 2

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimated reliability difference

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG, Wright TA (2015). “Cronbach's alpha reliability: Interval
estimation, hypothesis testing, and sample size planning.” *Journal of
Organizational Behavior*, **36**(1), 3–15. ISSN 08943796,
[doi:10.1002/job.1960](https://doi.org/10.1002/job.1960) . Bonett DG
(2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
ci.cronbach2(.05, .88, .76, 8, 8, 200, 250)
#>  Estimate     LL     UL
#>      0.12 0.0697 0.1733

# Should return:
# Estimate     LL     UL
#     0.12 0.0697 0.1733
 
```
