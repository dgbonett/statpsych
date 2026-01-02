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

Bonett DG, Wright TA (2015). “Cronbach's alpha reliability: Interval
estimation, hypothesis testing, and sample size planning.” *Journal of
Organizational Behavior*, **36**(1), 3–15. ISSN 08943796,
[doi:10.1002/job.1960](https://doi.org/10.1002/job.1960) .

## Examples

``` r
pi.score(.05, 24.5, 3.65, 40)
#>  Predicted df       LL       UL
#>       24.5 39 17.02546 31.97454

# Should return:
# Predicted  df       LL       UL
#      24.5  39 17.02546 31.97454
 
```
