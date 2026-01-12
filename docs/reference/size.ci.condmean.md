# Sample size for a conditional mean confidence interval

Computes the total sample size required to estimate a population
conditional mean of y at x = x\* in a fixed-x linear regression model
with desired confidence interval precision. The total sample size would
be allocated to the levels of the quantitative factor, and it might be
necessary to increase the total sample size to give the desired sample
size at each level of the fixed factor. Set the error variance planning
value to the largest value within a plausible range for a conservatively
large sample size.

For more details, see Section 1.24 of Bonett (2021, Volume 2)

## Usage

``` r
size.ci.condmean(alpha, evar, xvar, diff, w)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- evar:

  planning value of within group (error) variance

- xvar:

  variance of fixed predictor variable

- diff:

  difference between x\* and mean of x

- w:

  desired confidence interval width

## Value

Returns the required total sample size

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
size.ci.condmean(.05, 120, 125, 15, 5)
#>  Total sample size
#>                210

# Should return:
# Total sample size
#               210
 
```
