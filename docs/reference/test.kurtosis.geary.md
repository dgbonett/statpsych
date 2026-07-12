# Computes estimate and test of excess Geary kurtosis

Computes an estimate and test for kurtosis using a modfication of
Geary's measure of kurtosis. If the p-value is small (e.g., less than
.05) and excess kurtosis is positive, then the normality assumption can
be rejected due to leptokurtosis. If the p-value is small (e.g., less
than .05) and excess kurtosis is negative, then the normality assumption
can be rejected due to platykurtosis. The estimate and test of Geary's
kurtosis used here is based on a transformation of Geary's orginal
measure of kurtosis proposed by Bonett and Seier (2002). Geary's
kurtosis tends to be more sensitive to peakedness than Pearson's
kurtosis, and Pearson's kurtosis tend to be more sensitive to tail
weight than Geary's kurtosis. In the same way that it is informative to
assess centrality and variability using more than one measure, it is
also informative to assess kurtosis using both Pearson kurtosis and
Geary kurtosis. See (see
[test.kurtosis](https://dgbonett.github.io/statpsych/reference/test.kurtosis.md))
for a test of Pearson kurtosis.

## Usage

``` r
test.kurtosis.geary(y)
```

## Arguments

- y:

  vector of quantitative scores

## Value

Returns a 1-row matrix. The columns are:

- Kurtosis - estimate of transformed Geary kurtosis coefficient

- Excess - estimate of excess kurtosis (kurtosis - 3)

- z - z test statistic

- p - two-sided p-value

## References

Bonett DG, Seier E (2002). “A test of normality with high uniform
power.” *Computational Statistics & Data Analysis*, **40**, 435–445.
[doi:10.1016/S0167-9473(02)00074-9](https://doi.org/10.1016/S0167-9473%2802%2900074-9)
.

## Examples

``` r
y <- c(58, 58, 55, 52, 20, 65, 59, 49, 51, 81, 40, 62, 56, 49, 49, 50, 44, 53, 59)
test.kurtosis.geary(y)
#>  Kurtosis Excess     z     p
#>     5.157  2.157 2.793 0.005

# Should return:
# Kurtosis Excess     z     p
#    5.157  2.157 2.793 0.005

```
