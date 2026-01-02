# Confidence interval for a 2-group reliability difference

Computes a 100(1 - alpha)% confidence interval for a difference in
population reliabilities in a 2-group design. This function can be used
with any type of reliability coefficient (e.g., Cronbach alpha, McDonald
omega, intraclass reliability). The function requires a point estimate
and a 100(1 - alpha)% confidence interval for each reliability as input.

For more details, see Section 2.15 of Bonett (2021, Volume 4)

## Usage

``` r
ci.rel2(rel1, ll1, ul1, rel2, ll2, ul2)
```

## Arguments

- rel1:

  estimated reliability for group 1

- ll1:

  lower limit for group 1 reliability

- ul1:

  upper limit for group 1 reliability

- rel2:

  estimated reliability for group 2

- ll2:

  lower limit for group 2 reliability

- ul2:

  upper limit for group 2 reliability

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
ci.rel2(.4, .35, .47, .2, .1, .32)
#>  Estimate   LL     UL
#>       0.2 0.07 0.3221

# Should return:
# Estimate   LL     UL
#      0.2 0.07 0.3221
 
```
