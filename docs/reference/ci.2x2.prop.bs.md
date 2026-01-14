# Computes tests and confidence intervals of effects in a 2x2 between- subjects design for proportions

Computes adjusted Wald confidence intervals and tests for the AB
interaction effect, main effect of A, main effect of B, simple main
effects of A, and simple main effects of B in a 2x2 between-subjects
factorial design with a dichotomous response variable. The input vector
of frequency counts is f = \[ f11, f12, f21, f22 \], and the input
vector of sample sizes is n = \[ n11, n12, n21, n22 \] where the first
subscript represents the levels of Factor A and the second subscript
represents the levels of Factor B.

For more details, see Section 2.14 of Bonett (2021, Volume 3)

## Usage

``` r
ci.2x2.prop.bs(alpha, f, n)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- f:

  vector of frequency counts of participants who have the attribute

- n:

  vector of sample sizes

## Value

Returns a 7-row matrix (one row per effect). The columns are:

- Estimate - adjusted estimate of effect

- SE - standard error

- z - z test statistic

- p - two-sided p-value

- LL - lower limit of the adjusted Wald confidence interval

- UL - upper limit of the adjusted Wald confidence interval

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

Price RM, Bonett DG (2004). “An improved confidence interval for a
linear function of binomial proportions.” *Computational Statistics &
Data Analysis*, **45**(3), 449–456. ISSN 01679473,
[doi:10.1016/S0167-9473(03)00007-0](https://doi.org/10.1016/S0167-9473%2803%2900007-0)
.

## Examples

``` r
f <- c(15, 24, 28, 23)
n <- c(50, 50, 50, 50)
ci.2x2.prop.bs(.05, f, n)
#>             Estimate         SE       z       p          LL          UL
#> AB:      -0.27450980 0.13692496 -2.0048 0.04498 -0.54287780 -0.00614181
#> A:       -0.11764706 0.06846248 -1.7184 0.08572 -0.25183106  0.01653694
#> B:       -0.03921569 0.06846248 -0.5728 0.56678 -0.17339968  0.09496831
#> A at b1: -0.25000000 0.09402223 -2.6589 0.00784 -0.43428019 -0.06571981
#> A at b2:  0.01923077 0.09787658  0.1965 0.84422 -0.17260380  0.21106534
#> B at a1: -0.17307692 0.09432431 -1.8349 0.06652 -0.35794917  0.01179533
#> B at a2:  0.09615385 0.09758550  0.9853 0.32448 -0.09511021  0.28741790

# Should return:
#             Estimate         SE       z       p          LL          UL
# AB:      -0.27450980 0.13692496 -2.0048 0.04498 -0.54287780 -0.00614181
# A:       -0.11764706 0.06846248 -1.7184 0.08572 -0.25183106  0.01653694
# B:       -0.03921569 0.06846248 -0.5728 0.56678 -0.17339968  0.09496831
# A at b1: -0.25000000 0.09402223 -2.6589 0.00784 -0.43428019 -0.06571981
# A at b2:  0.01923077 0.09787658  0.1965 0.84423 -0.17260380  0.21106534
# B at a1: -0.17307692 0.09432431 -1.8349 0.06652 -0.35794917  0.01179533
# B at a2:  0.09615385 0.09758550  0.9853 0.32446 -0.09511021  0.28741790

```
