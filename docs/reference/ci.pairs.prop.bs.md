# Bonferroni confidence intervals for all pairwise proportion differences in a between-subjects design

Computes adjusted Wald confidence intervals for all pairwise differences
of population proportions in a between-subjects design using a
Bonferroni adjusted alpha level.

For more details, see Section 2.7 of Bonett (2021, Volume 3)

## Usage

``` r
ci.pairs.prop.bs(alpha, f, n)
```

## Arguments

- alpha:

  alpha level for simultaneous 1-alpha confidence

- f:

  vector of frequency counts of participants who have the attribute

- n:

  vector of sample sizes

## Value

Returns a matrix with the number of rows equal to the number of pairwise
comparisons. The columns are:

- Estimate - adjusted estimate of proportion difference

- SE - adjusted standard error

- z - z test statistic

- p - two-sided p-value

- LL - lower limit of the adjusted Wald confidence interval

- UL - upper limit of the adjusted Wald confidence interval

## References

Agresti A, Caffo B (2000). “Simple and effective confidence intervals
for proportions and differences of proportions result from adding two
successes and two failures.” *The American Statistician*, **54**(4),
280-288. ISSN 00031305,
[doi:10.2307/2685779](https://doi.org/10.2307/2685779) . Bonett DG
(2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
f <- c(111, 161, 132)
n <- c(200, 200, 200)
ci.pairs.prop.bs(.05, f, n)
#>        Estimate         SE       z       p          LL          UL
#>  1 2 -0.2475248 0.04482323 -5.5222 0.00000 -0.35483065 -0.14021885
#>  1 3 -0.1039604 0.04833562 -2.1508 0.03149 -0.21967489  0.01175409
#>  2 3  0.1435644 0.04358401  3.2940 0.00099  0.03922511  0.24790360

# Should return:
#        Estimate         SE       z       p          LL          UL
# 1 2  -0.2475248 0.04482323 -5.5222 0.00000 -0.35483065 -0.14021885
# 1 3  -0.1039604 0.04833562 -2.1508 0.03149 -0.21967489  0.01175409
# 2 3   0.1435644 0.04358401  3.2940 0.00099  0.03922511  0.24790360

```
