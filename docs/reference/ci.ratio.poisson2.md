# Confidence interval for a ratio of Poisson rates in a 2-group design

Computes a confidence interval for a ratio of population Poisson rates
in a 2-group design. The confidence interval is based on the binomial
method with an Agresti-Coull confidence interval. This function requires
the number of occurrences of a specific event (f) that were observed
over a specific period of time (t) within each group.

## Usage

``` r
ci.ratio.poisson2(alpha, f1, f2, t1, t2)
```

## Arguments

- alpha:

  alpha value for 1-alpha confidence

- f1:

  number of event occurrences for group 1

- f2:

  number of event occurrences for group 2

- t1:

  time period for group 1

- t2:

  time period for group 2

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimated ratio of Poisson rates

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## Details

The time periods do not need to be integers and can be expressed in any
unit of time such as seconds, hours, or months. The occurances are
assumed to be independent of one another and the unknown occurrence rate
is assumed to be constant over time within each group condition.

## References

Price RM, Bonett DG (2000). “Estimating the ratio of two Poisson rates.”
*Computational Statistics & Data Analysis*, **34**(3), 345–356.
[doi:10.1016/S0167-9473(99)00100-0](https://doi.org/10.1016/S0167-9473%2899%2900100-0)
.

## Examples

``` r
ci.ratio.poisson2(.05, 19, 5, 30, 40.5)
#>  Estimate       LL       UL
#>      5.13 1.939576 13.71481

# Should return:
# Estimate       LL       UL
#     5.13 1.939576 13.71481
 
```
