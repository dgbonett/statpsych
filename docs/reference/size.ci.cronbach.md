# Sample size for a Cronbach reliability confidence interval

Computes the sample size required to estimate a Cronbach reliability
with desired confidence interval precision. Set the reliability planning
value to the smallest value within a plausible range for a
conservatively large sample size.

## Usage

``` r
size.ci.cronbach(alpha, rel, r, w)
```

## Arguments

- alpha:

  alpha value for 1-alpha confidence

- rel:

  reliability planning value

- r:

  number of measurements (items, raters, forms)

- w:

  desired confidence interval width

## Value

Returns the required sample size

## Details

For more details, see Section 4.26 of Bonett (2021, Volume 1)

## References

Bonett DG, Wright TA (2015). “Cronbach's alpha reliability: Interval
estimation, hypothesis testing, and sample size planning.” *Journal of
Organizational Behavior*, **36**(1), 3–15. ISSN 08943796,
[doi:10.1002/job.1960](https://doi.org/10.1002/job.1960) .

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
size.ci.cronbach(.05, .75, 5, .15)
#>  Sample size
#>          109

# Should return:
# Sample size
#         109
 
```
