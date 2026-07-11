# Sample size to test an intraclass correlation

Computes the sample size required to test an intraclass correlation
(icc) with desired power in a one-way random effects ANOVA model.

## Usage

``` r
size.test.icc(alpha, pow, icc, r, h)
```

## Arguments

- alpha:

  alpha level for hypothesis test

- pow:

  desired power

- icc:

  icc planning value

- r:

  number of measurements (items, raters, forms)

- h:

  null hypothesis value of icc

## Value

Returns the required sample size

## References

Donner A, M. E (1987). “Sample size requirements for reliability
studies.” *Statistics in Medicine*, **6**, 441–448.
[doi:10.1002/sim.4780060404](https://doi.org/10.1002/sim.4780060404) .

## Examples

``` r
size.test.icc(.05, .90, .65, 4, .50)
#>  Sample size
#>          104

# Should return:
# Sample size
#         104
 
```
