# Test of a monotonic trend in means for an ordered between-subjects factor

Computes simultaneous confidence intervals for all adjacent pairwise
comparisons of population means using estimated group means, estimated
group standard deviations, and samples sizes as input. Equal variances
are not assumed. A Satterthwaite adjustment to the degrees of freedom is
used to improve the accuracy of the confidence intervals. If one or more
lower limits are greater than 0 and no upper limit is less than 0, then
conclude that the population means are monotonic decreasing. If one or
more upper limits are less than 0 and no lower limits are greater than
0, then conclude that the population means are monotonic increasing.
Reject the hypothesis of a monotonic trend if any lower limit is greater
than 0 and any upper limit is less than 0.

For more details, see Section 1.21 of Bonett (2021, Volume 2)

## Usage

``` r
test.mono.mean.bs(alpha, m, sd, n)
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

Returns a matrix with the number of rows equal to the number of adjacent
pairwise comparisons. The columns are:

- Estimate - estimated mean difference

- SE - standard error

- LL - one-sided lower limit of the confidence interval

- UL - one-sided upper limit of the confidence interval

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
m <- c(12.86, 24.57, 36.29, 53.21)
sd <- c(13.185, 12.995, 14.773, 15.145)
n <- c(20, 20, 20, 20)
test.mono.mean.bs(.05, m, sd, n)
#>      Estimate       SE        LL         UL
#>  1 2   -11.71 4.139530 -22.07803 -1.3419744
#>  2 3   -11.72 4.399497 -22.74731 -0.6926939
#>  3 4   -16.92 4.730817 -28.76921 -5.0707936

# Should return:
#     Estimate       SE        LL         UL
# 1 2   -11.71 4.139530 -22.07803 -1.3419744
# 2 3   -11.72 4.399497 -22.74731 -0.6926939
# 3 4   -16.92 4.730817 -28.76921 -5.0707936

```
