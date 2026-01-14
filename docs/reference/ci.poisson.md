# Confidence interval for a Poisson rate

Computes a confidence interval for a population Poisson rate. This
function requires the number of occurrences (f) of a specific event that
were observed over a specific period of time (t).

## Usage

``` r
ci.poisson(alpha, f, t)
```

## Arguments

- alpha:

  alpha value for 1-alpha confidence

- f:

  number of event occurrences

- t:

  time period

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimated Poisson rate

- SE - recovered standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## Details

The time period (t) does not need to be an integer and can be expressed
in any unit of time such as seconds, hours, or months. The occurances
are assumed to be independent of one another and the unknown occurance
rate is assumed to be constant over time.

## References

Hahn GJ, Meeker WQ (1991). *Statistical Intervals: A Guide for
Practitioners*. Wiley. ISBN 9780470316771,
[doi:10.1002/9780470316771](https://doi.org/10.1002/9780470316771) .

## Examples

``` r
ci.poisson(.05, 23, 5.25)
#>  Estimate        SE       LL      UL
#>  4.380952 0.9684952 2.777148 6.57358

# Should return:
# Estimate        SE       LL      UL
# 4.380952 0.9684952 2.777148 6.57358
 
```
