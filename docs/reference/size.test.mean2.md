# Sample size for a test of a 2-group mean difference

Computes the sample size in each group required to test a difference in
population means with desired power in a 2-group design. Set the
variance planning value to the largest value within a plausible range
for a conservatively large sample size. Set R =1 for equal sample sizes.
For unequal sample sizes, this function assumes approximately equal
population variances.

For more details, see Section 2.14 of Bonett (2021, Volume 1)

## Usage

``` r
size.test.mean2(alpha, pow, var, es, R)
```

## Arguments

- alpha:

  alpha level for hypothesis test

- pow:

  desired power

- var:

  planning value of average within-group variance

- es:

  planning value of mean difference

- R:

  n2/n1 ratio

## Value

Returns the required sample size for each group

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
size.test.mean2(.05, .95, 100, 10, 1) 
#>  n1 n2
#>  27 27

# Should return:
# n1  n2
# 27  27

size.test.mean2(.05, .95, 100, 10, 3) 
#>  n1 n2
#>  19 57

# Should return:
# n1  n2
# 19  57

size.test.mean2(.05, .95, 100, 10, .5) 
#>  n1 n2
#>  40 20

# Should return:
# n1  n2
# 40  20
 
```
