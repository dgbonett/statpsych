# Randomize a sample into groups

Randomly assigns a sample of participants into k groups.

## Usage

``` r
randomize(n)
```

## Arguments

- n:

  k x 1 vector of sample sizes

## Value

Returns a vector of randomly generated group assignments

## Examples

``` r
n <- c(10, 10, 5)
randomize(n)
#>  [1] 2 1 2 2 1 1 2 2 1 1 2 2 1 2 3 2 1 3 1 1 2 3 1 3 3

# Should return random numbers such as:
# [1] 2 3 2 1 1 2 3 3 2 1 2 1 3 1 1 2 3 1 1 2 2 1 1 2 2
 
```
