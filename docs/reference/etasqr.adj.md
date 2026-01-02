# Bias adjustment for an eta-squared estimate

Computes an approximate bias adjustment for eta-squared. This adjustment
can be applied to eta-squared, partial-eta squared, and generalized
eta-squared estimates.

For more details, see Section 3.7 of Bonett (2021, Volume 1)

## Usage

``` r
etasqr.adj(etasqr, dfeffect, dferror)
```

## Arguments

- etasqr:

  unadjusted eta-square estimate

- dfeffect:

  degrees of freedom for the effect

- dferror:

  error degrees of freedom

## Value

Returns a bias adjusted eta-squared estimate

## References

Bonett DG, Wright TA (2015). “Cronbach's alpha reliability: Interval
estimation, hypothesis testing, and sample size planning.” *Journal of
Organizational Behavior*, **36**(1), 3–15. ISSN 08943796,
[doi:10.1002/job.1960](https://doi.org/10.1002/job.1960) .

## Examples

``` r
etasqr.adj(.315, 2, 42)
#>  adj Eta-squared
#>           0.2824

# Should return:
# adj Eta-squared
#          0.2824
 
```
