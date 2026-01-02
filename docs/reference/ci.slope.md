# Confidence interval for a slope in a simple linear model

Computes a confidence interval for a population slope coefficient in a
simple linear model using the sample correlation, sample standard
deviation of the y scores (response variable), sample standard deviation
of the x scores (predictor variable), and sample size as input.

For more details, see Section 1.11 of Bonett (2021, Volume 2)

## Usage

``` r
ci.slope(alpha, cor, sdy, sdx, n)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- cor:

  estimated Pearson correlation

- sdy:

  estimated standard deviation of response variable

- sdx:

  estimated standard deviation of predictor variable

- n:

  sample size

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimated slope

- SE - standard error

- t - t test statistic

- df - degrees of freedom

- p - two-sided p-value

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
ci.slope(.05, .362, 25.1, 6.25, 85)
#>  Estimate        SE      t df       p        LL       UL
#>  1.453792 0.4109165 3.5379 83 0.00066 0.6364957 2.271088

# Should return:
#  Estimate        SE      t df       p        LL       UL
#  1.453792 0.4109165 3.5379 83 0.00066 0.6364957 2.271088
 
```
