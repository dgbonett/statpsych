# Generates random bivariate scores

Generates a random sample of y scores and x scores from a bivariate
normal distributions with specified population means, standard
deviations, and correlation. This function is useful for generating
hypothetical data for classroom demonstrations.

## Usage

``` r
random.yx(n, my, mx, sdy, sdx, cor, dec)
```

## Arguments

- n:

  sample size

- my:

  population mean of y scores

- mx:

  population mean of x scores

- sdy:

  population standard deviation of y scores

- sdx:

  population standard deviation of x scores

- cor:

  population correlation between x and y

- dec:

  number of decimal points

## Value

Returns n pairs of y and x scores

## Examples

``` r
random.yx(10, 50, 20, 4, 2, .5, 1)
#>       y    x
#> 1  46.0 21.2
#> 2  55.8 21.9
#> 3  41.0 17.8
#> 4  45.8 19.1
#> 5  50.3 18.8
#> 6  50.0 20.1
#> 7  49.2 19.9
#> 8  56.0 23.7
#> 9  47.6 20.7
#> 10 48.0 17.2

# Should return: 
#        y    x
#  1  50.3 21.6
#  2  52.0 21.6
#  3  53.0 22.7
#  4  46.9 21.3
#  5  56.3 23.8
#  6  50.4 20.3
#  7  44.6 19.9
#  8  49.9 18.3
#  9  49.4 18.5
# 10  42.3 20.2
 
```
