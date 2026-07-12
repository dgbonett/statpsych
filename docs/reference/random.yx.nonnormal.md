# Generates random bivariate nonnormal scores

Generates a random sample of y scores and x scores from a bivariate
nonnormal distribution with a specified population correlation and a
specified population mean, standard deviation, skewness, and excess
kurtosis for each variable. This function uses the mvrnonnorm function
in the semTools package.

For population excess kurtosis values greater than 0, the sample
kurtosis values tend to be smaller than the specified population values
because of their bias. The bias can be substanital for large excess
kurtosis values even in large samples.

## Usage

``` r
random.yx.nonnormal(n, my, mx, sdy, sdx, skewy, skewx, kury, kurx, cor, dec)
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

- skewy:

  population skewness of y scores

- skewx:

  population skewness of x scores

- kury:

  population excess kurtosis of y scores

- kurx:

  population excess kurtosis of x scores

- cor:

  population correlation between x and y

- dec:

  number of decimal points

## Value

Returns n pairs of y and x scores

## Examples

``` r
random.yx.nonnormal(10, 50, 20, 4, 2, .5, .75, 3, 4, .5, 1)
#>       y    x
#> 1  45.9 17.6
#> 2  47.8 19.6
#> 3  46.9 17.8
#> 4  49.3 16.5
#> 5  42.9 19.0
#> 6  49.4 20.2
#> 7  43.0 16.5
#> 8  44.7 17.2
#> 9  53.7 19.1
#> 10 47.1 18.9

# Should return:
#    y    x
# 1  47.3 18.8
# 2  52.0 20.3
# 3  51.3 22.3
# 4  50.9 21.3
# 5  55.2 22.0
# 6  53.7 20.1
# 7  46.4 20.0
# 8  54.9 24.7
# 9  48.7 20.7
# 10 44.7 18.6 
 
```
