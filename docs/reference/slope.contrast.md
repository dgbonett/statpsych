# Contrast coefficients for the slope of a quantitative factor

Computes the contrast coefficients that are needed to estimate the slope
of a line in a one-factor design with a quantitative factor.

For more details, see Section 1.11 of Bonett (2021, Volume 2)

## Usage

``` r
slope.contrast(x)
```

## Arguments

- x:

  vector of numeric factor levels

## Value

Returns the vector of contrast coefficients

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
x <- c(25, 50, 75, 100)
slope.contrast(x)
#>      Coefficient
#> [1,]      -0.012
#> [2,]      -0.004
#> [3,]       0.004
#> [4,]       0.012

# Should return:
#      Coefficient
# [1,]      -0.012
# [2,]      -0.004
# [3,]       0.004
# [4,]       0.012
 
```
