# Confidence interval for a Kendall tau-a correlation

Computes a Fisher confidence interval for a population Kendal tau-a
correlation. This function is useful for Kendall tau-a correlations less
than about .8. (Note: under bivariate normality, a Kendall tau-a
correlation of .8 corresponds to a Pearson correlation of about .95.) An
estimate of the Kendall tau-a correlation can be obtained from the (see
[cor.test](https://rdrr.io/r/stats/cor.test.html) function with method =
"kendall".

## Usage

``` r
ci.kendalltau(alpha, cor, n)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- cor:

  estimate of Kendall's tau-a correlation

- n:

  sample size

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimated Kendall tau-a correlation (from input)

- SE - standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG, Price RM (2002). “Statistical inference for a linear function
of medians: Confidence intervals, hypothesis testing, and sample size
requirements.” *Psychological Methods*, **7**(3), 370–383. ISSN
1939-1463,
[doi:10.1037/1082-989X.7.3.370](https://doi.org/10.1037/1082-989X.7.3.370)
.

## Examples

``` r
ci.kendalltau(.05, .734, 40)
#>  Estimate      SE     LL     UL
#>     0.734 0.05082 0.6178 0.8188

# Should return:
# Estimate      SE     LL     UL
#    0.734 0.05082 0.6178 0.8188

```
