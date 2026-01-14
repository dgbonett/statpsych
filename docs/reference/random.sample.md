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
#>  [1]    3  256  302  414  441  666  757  785  786  812  888  966 1150 1275 1348
#> [16] 2268 2360 2367 2396 2661 2774 2789 2915 2957 2960

# Should return random numbers such as:
#  [1]   37   94  134  186  212  408  485  697  722  781  998 1055 
# [13] 1182 1224 1273 1335 1452 1552 1783 1817 2149 2188 2437 2850 2936
 
```
