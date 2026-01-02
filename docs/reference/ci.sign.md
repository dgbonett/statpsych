# Confidence interval for the parameter of the one-sample sign test

Computes an adjusted Wald interval for the population proportion of
quantitative scores that are greater than the null hypothesis value of
the population median in a one-sample sign test. This proportion is a
measure of effect size that can be reported along with the sign test.

For more details, see Section 1.25 of Bonett (2021, Volume 1)

## Usage

``` r
ci.sign(alpha, y, h)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- y:

  vector of y scores

- h:

  null hypothesis value for population median

## Value

Returns a 1-row matrix. The columns are:

- Estimate - adjusted estimate of proportion

- SE - adjusted standard error

- LL - lower limit of adjusted Wald confidence interval

- UL - upper limit of adjusted Wald confidence interval

## References

Agresti A, Coull BA (1998). “Approximate is better than 'exact' for
interval estimation of binomial proportions.” *The American
Statistician*, **52**(2), 119–126. ISSN 0003-1305,
[doi:10.1080/00031305.1998.10480550](https://doi.org/10.1080/00031305.1998.10480550)
. Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
y <- c(30, 20, 15, 10, 10, 60, 20, 25, 20, 30, 10, 5, 50, 40, 20, 10,
        0, 20, 50)
ci.sign(.05, y, 9)
#>  Estimate        SE        LL        UL
#>  0.826087 0.0790342 0.6711828 0.9809911

# Should return:
# Estimate        SE        LL        UL
# 0.826087 0.0790342 0.6711828 0.9809911

```
