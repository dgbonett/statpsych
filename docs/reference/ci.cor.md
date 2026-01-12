# Confidence interval for a Pearson or partial correlation

Computes a Fisher confidence interval for a population Pearson
correlation or partial correlation with s control variables. Set s = 0
for a Pearson correlation. A bias adjustment is used to reduce the bias
of the Fisher transformed correlation. This function uses an estimated
correlation as input. Use the cor.test function for raw data input.

For more details, see Section 1.14 of Bonett (2021, Volume 2)

## Usage

``` r
ci.cor(alpha, cor, s, n)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- cor:

  estimated Pearson or partial correlation

- s:

  number of control variables

- n:

  sample size

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimated correlation (from input)

- SE - standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Snedecor GW, Cochran WG (1989). *Statistical Methods*, 8th edition. ISU
University Pres, Ames, Iowa.

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
ci.cor(.05, .60, 0, 150)
#>  Estimate      SE    LL     UL
#>       0.6 0.05243 0.485 0.6925

# Should return:
# Estimate      SE    LL     UL
#      0.6 0.05243 0.485 0.6925

ci.cor(.05, .70, 1, 135)
#>  Estimate      SE     LL     UL
#>       0.7 0.04406 0.6002 0.7763

# Should return:
# Estimate      SE     LL     UL
#      0.7 0.04406 0.6002 0.7763
 
```
