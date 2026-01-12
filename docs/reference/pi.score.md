# Prediction interval for one score

Computes a prediction interval for the response variable score of one
randomly selected member from the study population.

For more details, see Section 1.9 of Bonett (2021, Volume 1)

## Usage

``` r
pi.score(alpha, m, sd, n)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- m:

  estimated mean

- sd:

  estimated standard deviation

- n:

  sample size

## Value

Returns a 1-row matrix. The columns are:

- Predicted - predicted score

- df - degrees of freedom

- LL - lower limit of the prediction interval

- UL - upper limit of the prediction interval

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
pi.score(.05, 24.5, 3.65, 40)
#>  Predicted df       LL       UL
#>       24.5 39 17.02546 31.97454

# Should return:
# Predicted  df       LL       UL
#      24.5  39 17.02546 31.97454
 
```
