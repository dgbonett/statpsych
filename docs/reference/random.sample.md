# Generate a random sample

Generates a random sample of participant IDs without replacement.

## Usage

``` r
random.sample(popsize, samsize)
```

## Arguments

- popsize:

  study population size

- samsize:

  sample size

## Value

Returns a vector of randomly generated participant IDs

## Examples

``` r
random.sample(3000, 25)
#>  [1]   21   43  218  228  249  417  533  876  949  999 1051 1116 1372 1539 1546
#> [16] 1566 1615 1707 2093 2096 2201 2380 2395 2413 2807

# Should return random numbers such as:
#  [1]   37   94  134  186  212  408  485  697  722  781  998 1055 
# [13] 1182 1224 1273 1335 1452 1552 1783 1817 2149 2188 2437 2850 2936
 
```
