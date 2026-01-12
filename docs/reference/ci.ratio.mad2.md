# Confidence interval for a 2-group ratio of mean absolute deviations

Computes a confidence interval for a ratio of population MADs (mean
absolute deviation from median) in a 2-group design.

For more details, see Section 2.10 of Bonett (2021, Volume 1)

## Usage

``` r
ci.ratio.mad2(alpha, y1, y2)
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

- MAD1 - estimated MAD for group 1

- MAD2 - estimated MAD for group 2

- MAD1/MAD2 - estimate of MAD ratio

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG, Seier E (2003). “Confidence intervals for mean absolute
deviations.” *The American Statistician*, **57**(4), 233–236. ISSN
0003-1305,
[doi:10.1198/0003130032323](https://doi.org/10.1198/0003130032323) .

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
y1 <- c(32, 39, 26, 35, 43, 27, 40, 37, 34, 29)
y2 <- c(36, 44, 47, 42, 49, 39, 46, 31, 33, 48)
ci.ratio.mad2(.05, y1, y2)
#>      MAD1     MAD2 MAD1/MAD2        LL       UL
#>  5.111111 5.888889 0.8679245 0.4520879 1.666253

# Should return:
#     MAD1     MAD2  MAD1/MAD2        LL       UL
# 5.111111 5.888889  0.8679245 0.4520879 1.666253

```
