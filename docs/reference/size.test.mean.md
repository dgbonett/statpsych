# Sample size for a test of a mean

Computes the sample size required to test a single population mean with
desired power in a 1-group design. Set the variance planning value to
the largest value within a plausible range for a conservatively large
sample size.

For more details, see Section 1.29 of Bonett (2021, Volume 1)

## Usage

``` r
size.test.mean(alpha, pow, var, es)
```

## Arguments

- alpha:

  alpha level for hypothesis test

- pow:

  desired power

- var:

  planning value of response variable variance

- es:

  planning value of mean minus null hypothesis value

## Value

Returns the required sample size

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
size.test.mean(.05, .9, 8.2, 1.5)
#>  Sample size
#>           41

# Should return:
# Sample size
#          41
 
```
