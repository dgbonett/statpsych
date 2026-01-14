# Confidence interval for a Mann-Whitney parameter

Computes a distribution-free confidence interval for the Mann-Whitney
parameter (a "common language effect size"). In a 2-group experiment,
this parameter is the proportion of members in the population with
scores that would be higher under treatment 1 than treatment 2. In a
2-group nonexperiment where participants are sampled from two
subpopulations of sizes N1 and N2, the parameter is the proportion of
all N1 x N2 pairs in which a member from subpopulation 1 has a larger
score than a member from subpopulation 2.

For more details, see Section 2.12 of Bonett (2021, Volume 1)

## Usage

``` r
ci.mann(alpha, y1, y2)
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

- Estimate - estimated proportion

- SE - standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Sen PK (1967). “A note on asymptotically distribution-free confidence
bounds for P(X \< Y), based on two independent samples.” *The Indian
Journal of Statistics, Series A*, **29**(1), 95–102.

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
y1 <- c(9.4, 10.3, 58.3, 106.0, 31.0, 46.2, 12.0, 19.0, 135.0, 159.0)
y2 <- c(14.6, 5.1, 8.1, 22.7, 6.4, 4.4, 19.0, 3.2)
ci.mann(.05, y1, y2)
#>  Estimate        SE        LL UL
#>   0.86875 0.1222202 0.6292028  1

# Should return:
# Estimate        SE        LL UL
#  0.86875 0.1222202 0.6292028  1

```
