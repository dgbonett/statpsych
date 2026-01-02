# Sample size for a test of equal Pearson or partial correlation in a 2-group design

Computes the sample size required to test equality of two Pearson or
partial correlation with desired power in a 2-group design. Set s = 0
for a Pearson correlation. Set R = 1 for equal sample sizes.

## Usage

``` r
size.test.cor2(alpha, pow, cor1, cor2, s, R)
```

## Arguments

- alpha:

  alpha level for hypothesis test

- pow:

  desired power

- cor1:

  correlation planning value for group 1

- cor2:

  correlation planning value for group 2

- s:

  number of control variables

- R:

  n2/n1 ratio

## Value

Returns the required sample size for each group

## Examples

``` r
size.test.cor2(.05, .8, .4, .2, 0, 1)
#>   n1  n2
#>  325 325

# Should return:
#  n1  n2
# 325 325

size.test.cor2(.05, .8, .4, .2, 0, 2)
#>   n1  n2
#>  245 490

# Should return:
#  n1  n2
# 245 490
 
```
