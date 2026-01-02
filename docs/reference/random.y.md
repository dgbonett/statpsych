# Generate random sample of scores

Generates a random sample of scores from a normal distribution with a
specified population mean and standard deviation. This function is
useful for generating hypothetical data for classroom demonstrations.

## Usage

``` r
random.y(n, m, sd, min, max, dec)
```

## Arguments

- n:

  sample size

- m:

  population mean of scores

- sd:

  population standard deviation of scores

- min:

  minimum allowable value

- max:

  maximum allowable value

- dec:

  number of decimal points

## Value

Returns a vector of randomly generated scores.

## Examples

``` r
random.y(10, 3.6, 2.8, 1, 7, 0) 
#>  [1] 1 4 7 3 1 4 4 1 5 5

# Should return random numbers such as:
# [1] 2 7 7 1 6 3 1 3 2 1
 
```
