# Hypothesis test of equal proportions in a between-subjects design

Computes a Pearson chi-square test for equal population proportions for
a dichotomous response variable in a one-factor between-subjects design.

For more details, see Section 2.13 of Bonett (2021, Volume 3)

## Usage

``` r
test.prop.bs(f, n)
```

## Arguments

- f:

  vector of frequency counts of participants who have the attribute

- n:

  vector of sample sizes

## Value

Returns a 1-row matrix. The columns are:

- Chi-square - chi-square test statistic

- df - degrees of freedom

- p - p-value

## References

Fleiss JL, Paik MC (2003). *Statistical Methods for Rates and
Proportions*, 3rd edition. Wiley.

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
f <- c(111, 118, 132)
n <- c(200, 200, 200)
test.prop.bs (f, n)
#>  Chi-square df       p
#>      4.7706  2 0.09206

# Should return:
# Chi-square df       p
#     4.7706  2 0.09206

```
