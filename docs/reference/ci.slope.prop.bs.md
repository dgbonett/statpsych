# Confidence interval for a slope of a proportion in a single-factor experimental design with a quantitative between-subjects factor

Computes a test statistic and an adjusted Wald confidence interval for
the population slope of proportions in a one-factor experimental design
with a quantitative between-subjects factor.

For more details, see Section 4.4 of Bonett (2021, Volume 3)

## Usage

``` r
ci.slope.prop.bs(alpha, f, n, x)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- f:

  vector of frequency counts of participants who have the attribute

- n:

  vector of sample sizes

- x:

  vector of quantitative factor values

## Value

Returns a 1-row matrix. The columns are:

- Estimate - adjusted slope estimate

- SE - adjusted standard error

- z - z test statistic

- p - two-sided p-value

- LL - lower limit of the adjusted Wald confidence interval

- UL - upper limit of the adjusted Wald confidence interval

## References

Price RM, Bonett DG (2004). “An improved confidence interval for a
linear function of binomial proportions.” *Computational Statistics &
Data Analysis*, **45**(3), 449–456. ISSN 01679473,
[doi:10.1016/S0167-9473(03)00007-0](https://doi.org/10.1016/S0167-9473%2803%2900007-0)
. Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
f <- c(11, 15, 20, 27)
n <- c(60, 60, 60, 60)
x <- c(10, 20, 30, 40)
ci.slope.prop.bs(.05, f, n, x)
#>     Estimate          SE      z       p          LL         UL
#>  0.008688525 0.002566409 3.3855 0.00071 0.003658456 0.01371859

# Should return:
#    Estimate          SE      z       p          LL         UL
# 0.008688525 0.002566409 3.3855 0.00071 0.003658456 0.01371859

```
