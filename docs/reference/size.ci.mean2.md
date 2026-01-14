# Sample size for a 2-group mean difference confidence interval

Computes the sample size for each group required to estimate a
population mean difference with desired confidence interval precision in
a 2-group design. Set the variance planning value to the largest value
within a plausible range for a conservatively large sample size. Set R =
1 for equal sample sizes. For unequal sample sizes, this function
assumes approximately equal population variances.

For more details, see Section 2.13 of Bonett (2021, Volume 1)

## Usage

``` r
size.ci.mean2(alpha, var, w, R)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- var:

  planning value of average within-group variance

- w:

  desired confidence interval width

- R:

  n2/n1 ratio

## Value

Returns the required sample size for each group

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
size.ci.mean2(.01, 81, 12, 1)
#>  n1 n2
#>  32 32

# Should return:
# n1  n2
# 32  32

size.ci.mean2(.05, 4.0, 2.5, 2)
#>  n1 n2
#>  16 32

# Should return:
# n1  n2
# 16  32
 
```
