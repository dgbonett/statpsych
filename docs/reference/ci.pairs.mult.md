# Confidence intervals for pairwise proportion differences of a multinomial variable

Computes adjusted Wald confidence intervals for pairwise proportion
differences of a multinomial variable in a single sample. These adjusted
Wald confidence intervals use the same method that is used to compare
the two proportions in a paired-samples design.

For more details, see Section 1.12 of Bonett (2021, Volume 3)

## Usage

``` r
ci.pairs.mult(alpha, f)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- f:

  vector of multinomial frequency counts

## Value

Returns a matrix with the number of rows equal to the number of pairwise
comparisons. The columns are:

- Estimate - adjusted estimate of proportion difference

- SE - adjusted standard error

- LL - lower limit of the adjusted Wald confidence interval

- UL - upper limit of the adjusted Wald confidence interval

## References

Bonett DG, Price RM (2012). “Adjusted wald confidence interval for a
difference of binomial proportions based on paired data.” *Journal of
Educational and Behavioral Statistics*, **37**(4), 479–488. ISSN
1076-9986,
[doi:10.3102/1076998611411915](https://doi.org/10.3102/1076998611411915)
.

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
f <- c(87, 49, 31, 133)
ci.pairs.mult(.05, f)
#>         Estimate         SE           LL          UL
#>  1 2  0.12582781 0.03821865  0.050920628  0.20073500
#>  1 3  0.18543046 0.03466808  0.117482270  0.25337866
#>  1 4 -0.15231788 0.04855183 -0.247477718 -0.05715804
#>  2 3  0.05960265 0.02978792  0.001219398  0.11798590
#>  2 4 -0.27814570 0.04196760 -0.360400688 -0.19589070
#>  3 4 -0.33774834 0.03797851 -0.412184859 -0.26331183

# Should return:
#         Estimate         SE           LL          UL
#  1 2  0.12582781 0.03821865  0.050920628  0.20073500
#  1 3  0.18543046 0.03466808  0.117482270  0.25337866
#  1 4 -0.15231788 0.04855183 -0.247477718 -0.05715804
#  2 3  0.05960265 0.02978792  0.001219398  0.11798590
#  2 4 -0.27814570 0.04196760 -0.360400688 -0.19589070
#  3 4 -0.33774834 0.03797851 -0.412184859 -0.26331183

```
