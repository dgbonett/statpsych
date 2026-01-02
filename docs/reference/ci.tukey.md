# Tukey-Kramer confidence intervals for all pairwise mean differences in a between-subjects design

Computes heteroscedastic Tukey-Kramer (also known as Games-Howell)
confidence intervals for all pairwise comparisons of population means
using estimated means, estimated standard deviations, and samples sizes
as input. A Satterthwaite adjustment to the degrees of freedom is used
to improve the accuracy of the confidence intervals.

For more details, see Section 3.1 of Bonett (2021, Volume 1)

## Usage

``` r
ci.tukey(alpha, m, sd, n)
```

## Arguments

- alpha:

  alpha level for simultaneous 1-alpha confidence

- m:

  vector of estimated group means

- sd:

  vector of estimated group standard deviations

- n:

  vector of sample sizes

## Value

Returns a matrix with the number of rows equal to the number of pairwise
comparisons. The columns are:

- Estimate - estimated mean difference

- SE - standard error

- t - t test statistic

- df - degrees of freedom

- p - two-sided Tukey p-value

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Games PA, Howell JF (1976). “Pairwise multiple comparison procedures
with unequal N's and/or variances: A Monte Carlo study.” *Journal of
Educational Statistics*, **1**(2), 113. ISSN 03629791,
[doi:10.2307/1164979](https://doi.org/10.2307/1164979) . Bonett DG
(2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
m <- c(12.86, 17.57, 26.29, 30.21)
sd <- c(13.185, 12.995, 14.773, 15.145)
n <- c(20, 20, 20, 20)
ci.tukey(.05, m, sd, n)
#>      Estimate       SE       t    df       p        LL         UL
#>  1 2    -4.71 4.139530 -1.1378 37.99 0.66881 -15.83088  6.4108781
#>  1 3   -13.43 4.427673 -3.0332 37.52 0.02177 -25.33171 -1.5282918
#>  1 4   -17.35 4.490074 -3.8641 37.29 0.00233 -29.42285 -5.2771504
#>  2 3    -8.72 4.399497 -1.9820 37.39 0.21292 -20.54785  3.1078528
#>  2 4   -12.64 4.462292 -2.8326 37.14 0.03572 -24.64038 -0.6396178
#>  3 4    -3.92 4.730817 -0.8286 37.98 0.84056 -16.62952  8.7895242

# Should return:
#     Estimate       SE       t    df       p        LL         UL
# 1 2    -4.71 4.139530 -1.1378 37.99 0.66881 -15.83085  6.4108517
# 1 3   -13.43 4.427673 -3.0332 37.52 0.02177 -25.33172 -1.5282764
# 1 4   -17.35 4.490074 -3.8641 37.29 0.00233 -29.42281 -5.2771918
# 2 3    -8.72 4.399497 -1.9820 37.39 0.21291 -20.54783  3.1078269
# 2 4   -12.64 4.462292 -2.8326 37.14 0.03572 -24.64034 -0.6396589
# 3 4    -3.92 4.730817 -0.8286 37.98 0.84055 -16.62958  8.7895768

```
