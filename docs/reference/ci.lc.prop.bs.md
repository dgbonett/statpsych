# Confidence interval for a linear contrast of proportions in a between- subjects design

Computes an adjusted Wald confidence interval for a linear contrast of
population proportions in a between-subjects design.

For more details, see Section 2.8 of Bonett (2021, Volume 3)

## Usage

``` r
ci.lc.prop.bs(alpha, f, n, v)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- f:

  vector of frequency counts of participants who have the attribute

- n:

  vector of sample sizes

- v:

  vector of between-subjects contrast coefficients

## Value

Returns a 1-row matrix. The columns are:

- Estimate - adjusted estimate of proportion linear contrast

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
.

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
f <- c(26, 24, 38)
n <- c(60, 60, 60)
v1 <- c(-.5, -.5, 1)
ci.lc.prop.bs(.05/2, f, n, v1)
#>   Estimate         SE      z       p        LL        UL
#>  0.2119565 0.07602892 2.7878 0.00531 0.0415451 0.3823679

# Should return:
#  Estimate         SE      z       p        LL        UL
# 0.2119565 0.07602892 2.7878 0.00531 0.0415451 0.3823679

v2 <- c(1, -1, 0)
ci.lc.prop.bs(.05/2, f, n, v2)
#>    Estimate         SE      z       p         LL        UL
#>  0.03225806 0.08857951 0.3642 0.71571 -0.1662843 0.2308004

# Should return:
#   Estimate         SE       z       p         LL        UL
# 0.03225806 0.08857951 0.36417 0.71573 -0.1662843 0.2308004

```
