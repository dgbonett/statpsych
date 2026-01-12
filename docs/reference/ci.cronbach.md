# Confidence interval for a Cronbach reliability

Computes a confidence interval for a population Cronbach reliability.
The point estimate of Cronbach reliability assumes essentially
tau-equivalent measurements and the confidence interval assumes parallel
measurements.

For more details, see Section 4.19 of Bonett (2021, Volume 1)

## Usage

``` r
ci.cronbach(alpha, rel, r, n)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- rel:

  estimated Cronbach reliability

- r:

  number of measurements (items, raters, forms)

- n:

  sample size

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimated Cronbach reliability (from input)

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Feldt LS (1965). “The approximate sampling distribution of
Kuder-Richardson reliability coefficient twenty.” *Psychometrika*,
**30**(3), 357–370. ISSN 0033-3123,
[doi:10.1007/BF02289499](https://doi.org/10.1007/BF02289499) .

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
ci.cronbach(.05, .85, 7, 89)
#>  Estimate      SE     LL     UL
#>      0.85 0.02457 0.7971 0.8932

# Should return:
# Estimate       SE     LL     UL
#     0.85  0.02457 0.7971 0.8932
 
```
