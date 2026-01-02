# Confidence interval for a Spearman correlation

Computes a Fisher confidence interval for a population Spearman
correlation. Unlike the confidence interval for a Pearson correlation,
this function does not assume bivariate normality. Unlike the Pearson
correlation which describes a linear bivariate relation, the Spearman
correlation describes a monotonic bivariate relation. This function is
not appropriate for ordered categorical variables.

For more details, see Section 1.32 of Bonett (2021, Volume 2)

## Usage

``` r
ci.spear(alpha, y, x)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- y:

  vector of y scores

- x:

  vector of x scores (paired with y)

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimated Spearman correlation

- SE - standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG, Wright TA (2000). “Sample size requirements for estimating
Pearson, Kendall and Spearman correlations.” *Psychometrika*, **65**(1),
23–28. ISSN 0033-3123,
[doi:10.1007/BF02294183](https://doi.org/10.1007/BF02294183) . Bonett DG
(2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
y <- c(21, 4, 9, 12, 35, 18, 10, 22, 24, 1, 6, 8, 13, 16, 19)
x <- c(67, 28, 30, 28, 52, 40, 25, 37, 44, 10, 14, 20, 28, 40, 51)
ci.spear(.05, y, x)
#>  Estimate      SE     LL     UL
#>      0.87 0.08241 0.5841 0.9638

# Should return:
# Estimate      SE     LL     UL
#     0.87 0.08241 0.5841 0.9638
 
```
