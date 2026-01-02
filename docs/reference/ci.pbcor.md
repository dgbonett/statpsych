# Confidence intervals for point-biserial correlations

Computes confidence intervals for two types of population point-biserial
correlations. One type uses a weighted average of the group variances
and is appropriate for nonexperimental designs with simple random
sampling (but not stratified random sampling). The other type uses an
unweighted average of the group variances and is appropriate for
experimental designs. Equality of variances is not assumed for either
type.

For more details, see Section 1.18 of Bonett (2021, Volume 2)

## Usage

``` r
ci.pbcor(alpha, m1, m2, sd1, sd2, n1, n2)
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

## Value

Returns a 2-row matrix. The columns are:

- Estimate - estimated point-biserial correlation

- SE - standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG (2020). “Point-biserial correlation: Interval estimation,
hypothesis testing, meta-analysis, and sample size determination.”
*British Journal of Mathematical and Statistical Psychology*,
**73**(S1), 113–144. ISSN 0007-1102,
[doi:10.1111/bmsp.12189](https://doi.org/10.1111/bmsp.12189) . Bonett DG
(2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
ci.pbcor(.05, 28.32, 21.48, 3.81, 3.09, 40, 40)
#>             Estimate      SE     LL     UL
#> Weighted:     0.7066 0.04891 0.5885 0.7854
#> Unweighted:   0.7021 0.05019 0.5808 0.7829

# Should return:
#             Estimate      SE     LL     UL
# Weighted:     0.7066 0.04891 0.5885 0.7854
# Unweighted:   0.7021 0.05019 0.5808 0.7829
 
```
