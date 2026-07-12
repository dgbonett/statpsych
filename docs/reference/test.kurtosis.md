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
y <- c(58, 58, 55, 52, 20, 65, 59, 49, 51, 81, 40, 62, 56, 49, 49, 50, 44, 53, 59)
test.kurtosis(y)
#>  Kurtosis Excess     p
#>     5.508  2.508 0.017

# Should return:
# Kurtosis Excess      p
#   5.5081 2.5081 0.0164

```
