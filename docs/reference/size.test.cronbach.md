# Sample size to test a Cronbach reliability

Computes the sample size required to test a Cronbach reliability with
desired power.

## Usage

``` r
size.test.cronbach(alpha, pow, rel, r, h)
```

## Arguments

- alpha:

  alpha level for hypothesis test

- pow:

  desired power

- rel:

  reliability planning value

- r:

  number of measurements (items, raters, forms)

- h:

  null hypothesis value of reliability

## Value

Returns the required sample size

## References

Bonett DG, Wright TA (2015). “Cronbach's alpha reliability: Interval
estimation, hypothesis testing, and sample size planning.” *Journal of
Organizational Behavior*, **36**(1), 3–15. ISSN 08943796,
[doi:10.1002/job.1960](https://doi.org/10.1002/job.1960) .

## Examples

``` r
size.test.cronbach(.05, .85, .80, 5, .7)
#>  Sample size
#>          139

# Should return:
# Sample size
#         139
 
```
