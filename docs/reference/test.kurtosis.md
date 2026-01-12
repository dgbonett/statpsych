# Computes p-value for test of excess kurtosis

Computes a Monte Carlo p-value (250,000 replications) for the null
hypothesis that the sample data come from a normal distribution. If the
p-value is small (e.g., less than .05) and excess kurtosis is positive,
then the normality assumption can be rejected due to leptokurtosis. If
the p-value is small (e.g., less than .05) and excess kurtosis is
negative, then the normality assumption can be rejected due to
platykurtosis.

For more details, see Section 1.23 of Bonett (2021, Volume 1)

## Usage

``` r
test.kurtosis(y)
```

## Arguments

- y:

  vector of quantitative scores

## Value

Returns a 1-row matrix. The columns are:

- Kurtosis - estimate of kurtosis coefficient

- Excess - estimate of excess kurtosis (kurtosis - 3)

- p - Monte Carlo two-sided p-value

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
y <- c(30, 20, 15, 10, 10, 60, 20, 25, 20, 30, 10, 5, 50, 40, 95)
test.kurtosis(y)
#>  Kurtosis Excess    p
#>     4.815  1.815 0.04

# Should return:
# Kurtosis  Excess     p
#   4.8149  1.8149 0.038 

```
