# Sample size for a slope confidence interval

Computes the total sample size required to estimate a population slope
with desired confidence interval precision in a between-subjects design
with a quantitative factor. In an experimental design, the total sample
size would be allocated to the levels of the quantitative factor and it
might be necessary to increase the total sample size to achieve equal
sample sizes. Set the error variance planning value to the largest value
within a plausible range for a conservatively large sample size.

For more details, see Section 1.24 of Bonett (2021, Volume 2)

## Usage

``` r
size.ci.slope(alpha, evar, x, w)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- evar:

  planning value of within-group (error) variance

- x:

  vector of x values of the quantitative factor

- w:

  desired confidence interval width

## Value

Returns the required total sample size

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
x <- c(5, 7, 9)
size.ci.slope(.05, 300, x, 5)
#>  Total sample size
#>                 73

# Should return:
# Total sample size
#                73
 
```
