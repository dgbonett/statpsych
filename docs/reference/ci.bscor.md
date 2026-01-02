# Confidence interval for a biserial correlation

Computes a confidence interval for a population biserial correlation. A
biserial correlation can be used when one variable is quantitative and
the other variable has been artificially dichotomized to create two
groups. The biserial correlation estimates the correlation between the
observed quantitative variable and the unobserved quantitative variable
that has been measured on a dichotomous scale.

## Usage

``` r
ci.bscor(alpha, m1, m2, sd1, sd2, n1, n2)
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

Returns a 1-row matrix. The columns are:

- Estimate - estimated biserial correlation

- SE - standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## Details

This function computes a point-biserial correlation and its standard
error as a function of a standardized mean difference with a weighted
variance standardizer. Then the point-biserial estimate is transformed
into a biserial correlation using the traditional adjustment. The
adjustment is also applied to the point-biserial standard error to
obtain the standard error for the biserial correlation.

The biserial correlation assumes that the observed quantitative variable
and the unobserved quantitative variable have a bivariate normal
distribution. Bivariate normality is a crucial assumption underlying the
transformation of a point-biserial correlation to a biserial
correlation. Bivariate normality also implies equal variances of the
observed quantitative variable at each level of the dichotomized
variable, and this assumption is made in the computation of the standard
error.

## References

Bonett DG (2020). “Point-biserial correlation: Interval estimation,
hypothesis testing, meta-analysis, and sample size determination.”
*British Journal of Mathematical and Statistical Psychology*,
**73**(S1), 113–144. ISSN 0007-1102,
[doi:10.1111/bmsp.12189](https://doi.org/10.1111/bmsp.12189) .

## Examples

``` r
ci.bscor(.05, 28.32, 21.48, 3.81, 3.09, 40, 40)
#>  Estimate     SE     LL     UL
#>    0.8856 0.0613 0.7376 0.9844

# Should return:
# Estimate     SE      LL     UL
#   0.8856 0.0613  0.7376 0.9844
 
```
