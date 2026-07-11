# Generates random bivariate normal scores

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
#> 1  51.2 18.6
#> 2  50.4 19.4
#> 3  49.0 16.2
#> 4  48.9 17.8
#> 5  45.7 18.4
#> 6  48.7 20.4
#> 7  51.3 21.3
#> 8  48.0 19.2
#> 9  41.9 15.0
#> 10 43.5 18.2

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
