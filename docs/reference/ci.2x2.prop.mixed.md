# Computes tests and confidence intervals of effects in a 2x2 mixed design for proportions

Computes adjusted Wald confidence intervals and tests for the AB
interaction effect, main effect of A, main effect of B, simple main
effects of A, and simple main effects of B in a 2x2 mixed factorial
design with a dichotomous response variable where Factor A is a
within-subjects factor and Factor B is a between-subjects factor. The
4x1 vector of frequency counts for Factor A within each group is f00,
f01, f10, f11 where fij is the number of participants with a response of
i = 0 or 1 at level 1 of Factor A and a response of j = 0 or 1 at level
2 of Factor A.

For more details, see Section 2.2 of Bonett (2021, Volume 3)

## Usage

``` r
ci.2x2.prop.mixed(alpha, group1, group2)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- group1:

  vector of frequency counts from 2x2 contingency table in group 1

- group2:

  vector of frequency counts from 2x2 contingency table in group 2

## Value

Returns a 7-row matrix (one row per effect). The columns are:

- Estimate - adjusted estimate of effect

- SE - standard error of estimate

- z - z test statistic

- p - two-sided p-value

- LL - lower limit of the adjusted Wald confidence interval

- UL - upper limit of the adjusted Wald confidence interval

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
group1 <- c(125, 14, 10, 254)
group2 <- c(100, 16, 9, 275)
ci.2x2.prop.mixed (.05, group1, group2)
#>              Estimate          SE       z       p          LL           UL
#> AB:       0.007555369 0.017716073  0.4265 0.66974 -0.02716750  0.042278234
#> A:       -0.013678675 0.008858036 -1.5442 0.12254 -0.03104011  0.003682758
#> B:       -0.058393219 0.023032656 -2.5352 0.01124 -0.10353640 -0.013250043
#> A at b1: -0.009876543 0.012580603 -0.7851 0.43239 -0.03453407  0.014780985
#> A at b2: -0.017412935 0.012896543 -1.3502 0.17695 -0.04268969  0.007863824
#> B at a1: -0.054634236 0.032737738 -1.6688 0.09516 -0.11879902  0.009530550
#> B at a2: -0.062170628 0.032328556 -1.9231 0.05447 -0.12553343  0.001192177

# Should return:
#              Estimate          SE       z       p          LL           UL
# AB:       0.007555369 0.017716073  0.4265 0.66976 -0.02716750  0.042278234
# A:       -0.013678675 0.008858036 -1.5442 0.12254 -0.03104011  0.003682758
# B:       -0.058393219 0.023032656 -2.5352 0.01124 -0.10353640 -0.013250043
# A at b1: -0.009876543 0.012580603 -0.7851 0.43242 -0.03453407  0.014780985
# A at b2: -0.017412935 0.012896543 -1.3502 0.17695 -0.04268969  0.007863824
# B at a1: -0.054634236 0.032737738 -1.6688 0.09515 -0.11879902  0.009530550
# B at a2: -0.062170628 0.032328556 -1.9231 0.05447 -0.12553343  0.001192177

```
