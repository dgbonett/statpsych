# Between-subjects F statistic and eta-squared from summary information

Computes the F statistic, p-value, eta-squared, and adjusted eta-squared
for the main effect in a one-way between-subjects ANOVA using the
estimated group means, estimated group standard deviations, and group
sample sizes.

For more details, see Section 3.7 of Bonett (2021, Volume 1)

## Usage

``` r
test.anova.bs(m, sd, n)
```

## Arguments

- m:

  vector of estimated group means

- sd:

  vector of estimated group standard deviations

- n:

  vector of group sample sizes

## Value

Returns a 1-row matrix. The columns are:

- F - F statistic for test of null hypothesis

- dfA - degrees of freedom for between-subjects factor

- dfE - error degrees of freedom

- p - p-value

- Eta-squared - estimate of eta-squared

- adj Eta-squared - a bias adjusted estimate of eta-squared

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
m <- c(12.4, 8.6, 10.5)
sd <- c(3.84, 3.12, 3.48)
n <- c(20, 20, 20)
test.anova.bs(m, sd, n)
#>       F dfA dfE       p Eta-squared adj Eta-squared
#>  5.9196   2  57 0.00461       0.172          0.1429

#  Should return:
#      F dfA  dfE       p Eta-squared  adj Eta-squared
# 5.9196   2   57 0.00461       0.172           0.1429
 
```
