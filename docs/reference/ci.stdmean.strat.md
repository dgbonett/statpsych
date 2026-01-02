# Confidence intervals for a 2-group standardized mean difference with stratified sampling

Computes confidence intervals for a population standardized mean
difference in a 2-group nonexperimental design with stratified random
sampling (a random sample of a specified size from each subpopulation)
using a square root weighted variance standardizer or single group
standard deviation standardizer. Equality of variances is not assumed.

## Usage

``` r
ci.stdmean.strat(alpha, m1, m2, sd1, sd2, n1, n2, p1)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- m1:

  estimated mean for group 1

- m2:

  estimated mean for group 2

- sd1:

  estimated standard deviation for group 1

- sd2:

  estimated standard deviation for group 2

- n1:

  sample size for group 1

- n2:

  sample size for group 2

- p1:

  proportion of total population in subpopulation 1

## Value

Returns a 3-row matrix. The columns are:

- Estimate - estimated standardized mean difference

- adj Estimate - bias adjusted standardized mean difference estimate

- SE - standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG (2020). “Point-biserial correlation: Interval estimation,
hypothesis testing, meta-analysis, and sample size determination.”
*British Journal of Mathematical and Statistical Psychology*,
**73**(S1), 113–144. ISSN 0007-1102,
[doi:10.1111/bmsp.12189](https://doi.org/10.1111/bmsp.12189) .

## Examples

``` r
ci.stdmean.strat(.05, 33.2, 30.8, 10.5, 11.2, 200, 200, .533)
#>                        Estimate adj Estimate      SE     LL     UL
#> Weighted standardizer:   0.2216       0.2211 0.10052 0.0245 0.4186
#> Group 1 standardizer:    0.2286       0.2277 0.10428 0.0242 0.4330
#> Group 2 standardizer:    0.2143       0.2277 0.09776 0.0227 0.4059

# Should return:
#                        Estimate adj Estimate      SE     LL     UL
# Weighted standardizer:   0.2216       0.2211 0.10052 0.0245 0.4186
# Group 1 standardizer:    0.2286       0.2277 0.10428 0.0242 0.4330
# Group 2 standardizer:    0.2143       0.2277 0.09776 0.0227 0.4059

```
