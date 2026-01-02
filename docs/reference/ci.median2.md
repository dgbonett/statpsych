# Confidence interval for a 2-group median difference

Computes a distribution-free confidence interval for a difference of
population medians in a 2-group design. Tied scores within each group
are assumed to be rare.

For more details, see Section 2.12 of Bonett (2021, Volume 1)

## Usage

``` r
ci.median2(alpha, y1, y2)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- y1:

  vector of scores for group 1

- y2:

  vector of scores for group 2

## Value

Returns a 1-row matrix. The columns are:

- Median1 - estimated median for group 1

- Median2 - estimated median for group 2

- Median1-Median2 - estimated difference of medians

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
y1 = c(70, 394, 43, 95, 62, 128, 2, 203, 81, 436, 85, 35, 156, 1, 3, 27, 63, 181, 184, 18)
y2 = c(102, 120, 78, 78, 417, 124, 86, 171, 176, 129, 230,194, 100, 157, 306, 411, 164, 103, 193, 312)
ci.median2(.05, y1, y2)
#>  Median1 Median2 Median1-Median2       SE        LL        UL
#>     75.5   160.5             -85 37.11502 -157.7441 -12.25589

# Should return:
#  Median1 Median2 Median1-Median2       SE        LL        UL
#     75.5   160.5             -85 37.11502 -157.7441 -12.25589

```
