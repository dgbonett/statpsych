# Confidence intervals for diversity indices

Computes estimates, standard errors, and approximate confidence
intervals for the Berger-Parker, Simpson, and Shannon diversity indices.
For the Shannon index, the value 1/r is added to each frequency count
where r is the number of categories. These indices have a range of 0 to
1 where 0 indicates no diversity and 1 indicates maximum diversity.

For more details, see Section 1.13 of Bonett (2021, Volume 3)

## Usage

``` r
ci.diversity(alpha, f)
```

## Arguments

- alpha:

  alpha level for 1 - alpha confidence

- f:

  vector of multinomial frequency counts

## Value

Returns a 3-row matrix. The columns are:

- Estimate - estimate of diversity index

- SE - standard error of estimate

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
f = c(847, 320, 57, 274, 36)
ci.diversity(.05, f)
#>         Estimate      SE     LL     UL
#> Berger    0.5598 0.01587 0.5287 0.5909
#> Simpson   0.7722 0.01229 0.7481 0.7963
#> Shannon   0.7292 0.01224 0.7052 0.7532

# Should return:
#         Estimate      SE     LL     UL
# Berger    0.5598 0.01587 0.5287 0.5909
# Simpson   0.7722 0.01229 0.7481 0.7963
# Shannon   0.7292 0.01224 0.7052 0.7532
 
```
