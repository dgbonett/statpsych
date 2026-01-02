# Confidence interval for a linear contrast of medians in a between-subjects design

Computes a distribution-free confidence interval for a linear contrast
of medians in a between-subjects design using estimated medians and
their standard errors. The sample median and standard error for each
group can be computed using the
[ci.median](https://dgbonett.github.io/statpsych/reference/ci.median.md)
function.

For more details, see Section 3.21 of Bonett (2021, Volume 1)

## Usage

``` r
ci.lc.median.bs(alpha, m, se, v)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- m:

  vector of estimated group medians

- se:

  vector of estimated group standard errors

- v:

  vector of between-subjects contrast coefficients

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimated linear contrast of medians

- SE - standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG, Price RM (2002). “Statistical inference for a linear function
of medians: Confidence intervals, hypothesis testing, and sample size
requirements.” *Psychological Methods*, **7**(3), 370–383. ISSN
1939-1463,
[doi:10.1037/1082-989X.7.3.370](https://doi.org/10.1037/1082-989X.7.3.370)
. Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
m <- c(46.13, 29.19, 30.32, 49.15)
se <- c(6.361, 5.892, 4.887, 6.103)
v <- c(1, -1, -1, 1)
ci.lc.median.bs(.05, m, se, v)
#>  Estimate       SE       LL       UL
#>     35.77 11.67507 12.88727 58.65273

# Should return:
# Estimate       SE       LL       UL
#    35.77 11.67507 12.88727 58.65273

```
