# Confidence interval for the slope of means in a one-factor experimental design with a quantitative between-subjects factor

Computes a test statistic and confidence interval for the slope of means
in a one-factor experimental design with a quantitative between-subjects
factor. This function computes both the unequal variance and equal
variance confidence intervals and test statistics. A Satterthwaite
adjustment to the degrees of freedom is used with the unequal variance
method.

For more details, see Section 1.11 of Bonett (2021, Volume 2)

## Usage

``` r
ci.slope.mean.bs(alpha, m, sd, n, x)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- m:

  vector of sample means

- sd:

  vector of sample standard deviations

- n:

  vector of sample sizes

- x:

  vector of quantiative factor values

## Value

Returns a 2-row matrix. The columns are:

- Estimate - estimated slope

- SE - standard error

- t - t test statistic

- df - degrees of freedom

- p - two-sided p-value

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
m <- c(33.5, 37.9, 38.0, 44.1)
sd <- c(3.84, 3.84, 3.65, 4.98)
n <- c(10,10,10,10)
x <- c(5, 10, 20, 30)
ci.slope.mean.bs(.05, m, sd, n, x)
#>                               Estimate         SE      t    df     p        LL
#> Equal Variances Assumed:     0.3664407 0.06770529 5.4123 36.00 0e+00 0.2291280
#> Equal Variances Not Assumed: 0.3664407 0.07336289 4.9949 18.66 8e-05 0.2127008
#>                                     UL
#> Equal Variances Assumed:     0.5037534
#> Equal Variances Not Assumed: 0.5201806

# Should return:
#                               Estimate         SE      t    df
# Equal Variances Assumed:     0.3664407 0.06770529 5.4123 36.00
# Equal Variances Not Assumed: 0.3664407 0.07336289 4.9949 18.66
#                                  p        LL        UL
# Equal Variances Assumed:     4e-06 0.2291280 0.5037534
# Equal Variances Not Assumed: 8e-05 0.2126998 0.5201815

```
