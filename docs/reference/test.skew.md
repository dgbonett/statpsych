# Computes p-value for test of skewness

Computes a Monte Carlo p-value (250,000 replications) for the null
hypothesis that the sample data come from a normal distribution. If the
p-value is small (e.g., less than .05) and the skewness estimate is
positive, then the normality assumption can be rejected due to positive
skewness. If the p-value is small (e.g., less than .05) and the skewness
estimate is negative, then the normality assumption can be rejected due
to negative skewness.

For more details, see Section 1.23 of Bonett (2021, Volume 1)

## Usage

``` r
test.skew(y)
```

## Arguments

- y:

  vector of quantitative scores

## Value

Returns a 1-row matrix. The columns are:

- Skewness - estimate of skewness coefficient

- p - Monte Carlo two-sided p-value

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
y <- c(30, 20, 15, 10, 10, 60, 20, 25, 20, 30, 10, 5, 50, 40, 95)
test.skew(y)
#>  Skewness     p
#>      1.52 0.007

# Should return:
# Skewness     p
#   1.5201 0.007

```
