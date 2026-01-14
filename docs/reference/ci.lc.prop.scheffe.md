# Scheffe confidence interval for a linear contrast of proportions in a between-subjects design

Computes an adjusted Wald confidence interval for a linear contrast of
population proportions in a between-subjects design using a Scheffe
critical value. A Scheffe p-value is computed for the test statistic.
This function is useful in exploratory studies where the linear contrast
of proportions was not planned but was suggested by the pattern of
sample proportions. Use the
[ci.lc.prop.bs](https://dgbonett.github.io/statpsych/reference/ci.lc.prop.bs.md)
function with a Bonferroni adjusted alpha value to compute simultaneous
confidence intervals for two or more planned linear contrasts of
proportions.

For more details, see Section 2.9 of Bonett (2021, Volume 3)

## Usage

``` r
ci.lc.prop.scheffe(alpha, f, n, v)
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

- p - two-sided Scheffe p-value

- LL - lower limit of the Scheffe confidence interval

- UL - upper limit of the Scheffe confidence interval

## References

Price RM, Bonett DG (2004). “An improved confidence interval for a
linear function of binomial proportions.” *Computational Statistics &
Data Analysis*, **45**(3), 449–456. ISSN 01679473,
[doi:10.1016/S0167-9473(03)00007-0](https://doi.org/10.1016/S0167-9473%2803%2900007-0)
.

Marascuilo LA, McSweeney M (1977). *Nonparametric and Distribution-Free
Methods for the Social Sciences*. Brooks/Cole.

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
f <- c(26, 24, 38)
n <- c(60, 60, 60)
v <- c(-.5, -.5, 1)
ci.lc.prop.scheffe(.05, f, n, v)
#>   Estimate         SE      z       p         LL        UL
#>  0.2119565 0.07602892 2.7878 0.02053 0.02585698 0.3980561

# Should return:
#  Estimate         SE      z       p         LL        UL
# 0.2119565 0.07602892 2.7878 0.02053 0.02585698 0.3980561

```
