# Prediction interval for difference of scores in a 2-level within-subjects experiment

For a 2-level within-subjects experiment, this function computes a
prediction interval for how the response variable score for one randomly
selected person from the study population would differ under the two
treatment conditions.

For more details, see Section 4.5 of Bonett (2021, Volume 1)

## Usage

``` r
pi.score.ps(alpha, m1, m2, sd1, sd2, cor, n)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- m1:

  estimated mean for measurement 1

- m2:

  estimated mean for measurement 2

- sd1:

  estimated standard deviation for measurement 1

- sd2:

  estimated standard deviation for measurement 2

- cor:

  estimated correlation of paired scores

- n:

  sample size

## Value

Returns a 1-row matrix. The columns are:

- Predicted - predicted difference in scores

- df - degrees of freedom

- LL - lower limit of the prediction interval

- UL - upper limit of the prediction interval

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
pi.score.ps(.05, 265.1, 208.6, 23.51, 19.94, .814, 30)
#>  Predicted df       LL       UL
#>       56.5 29 28.05936 84.94064

# Should return:
# Predicted df       LL       UL
#      56.5 29 28.05936 84.94064
 
```
