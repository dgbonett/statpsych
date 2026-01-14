# Sample size for a intraclass correlation confidence interval

Computes the sample size required to estimate an intraclass correlation
with desired confidence interval precision. This type of intraclass
correlation can be used to describe the reliability of a single
measurement (e.g., a single rater or a single form of a test). This
intraclass correlation assumes a two-factor (subject x measurement)
model. Use the size.ci.cronbach function to determine the sample size
required to estimate the reliability of a sum or average of r
measurements (e.g., items) with desired confidence interval precision.

## Usage

``` r
size.ci.icc(alpha, icc, r, w)
```

## Arguments

- alpha:

  alpha value for 1-alpha confidence

- icc:

  intraclass correlation planning value

- r:

  number of measurements (raters, forms, occasions)

- w:

  desired confidence interval width

## Value

Returns the required sample size

## Details

Specifying an intraclass correlation planning value for a conservatively
large sample size requirement is not straightforward with an intraclass
correlation. For r = 2, the sample size requirement is largest for an
intraclass correlation of zero. But for r = 3, 4, 5, 10, 20, and 40 an
intraclass correlation planning value of about .26, .33, .38, .43, .44,
and .45, respectively, maximizes the sample size requirement.

## References

Bonett DG (2002). “Sample size requirements for estimating intraclass
correlations with desired precision.” *Statistics in Medicine*, **21**,
1331–1335. [doi:10.1002/sim.1108](https://doi.org/10.1002/sim.1108) .

## Examples

``` r
size.ci.icc(.05, .70, 3, .2)
#>  Sample size
#>           68

# Should return:
# Sample size
#          68
 
```
