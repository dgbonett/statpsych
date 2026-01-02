# Sample size for a test of a slope

Computes the total sample size required to test a population slope with
desired power in a between-subjects design with a quantitative factor.
In an experimental design, the total sample size would be allocated to
the levels of the quantitative factor and it might be necessary to use a
larger total sample size to achieve equal sample sizes. Set the error
variance planning value to the largest value within a plausible range
for a conservatively large sample size.

For more details, see Section 1.25 of Bonett (2021, Volume 2)

## Usage

``` r
size.test.slope(alpha, pow, evar, x, slope, h)
```

## Arguments

- alpha:

  alpha level for hypothesis test

- pow:

  desired power

- evar:

  planning value of within-group (error) variance

- x:

  vector of x values of the quantitative factor

- slope:

  planning value of slope

- h:

  null hypothesis value of slope

## Value

Returns the required total sample size

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
x <- c(2, 5, 8)
size.test.slope(.05, .9, 31.1, x, .75, 0)
#>  Total sample size
#>                100

# Should return:
# Total sample size
#               100
 
```
