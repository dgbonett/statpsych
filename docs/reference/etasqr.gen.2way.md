# Generalized eta-squared estimates in a two-factor design

Computes generalized eta-square estimates in a two-factor design where
one or both factors are classification factors. If both factors are
treatment factors, then partial eta-square estimates are typically
recommended. The eta-squared estimates from this function can be used in
the
[etasqr.adj](https://dgbonett.github.io/statpsych/reference/etasqr.adj.md)
function to obtain bias adjusted estimates.

For more details, see Section 3.12 of Bonett (2021, Volume 1)

## Usage

``` r
etasqr.gen.2way(SSa, SSb, SSab, SSe)
```

## Arguments

- SSa:

  sum of squares for factor A

- SSb:

  sum of squares for factor B

- SSab:

  sum of squares for A x B interaction

- SSe:

  error (within) sum of squares

## Value

Returns a 3-row matrix. The columns are:

- A - estimate of eta-squared for factor A

- B - estimate of eta-squared for factor B

- AB - estimate of eta-squared for A x B interaction

## References

Bonett DG, Wright TA (2015). “Cronbach's alpha reliability: Interval
estimation, hypothesis testing, and sample size planning.” *Journal of
Organizational Behavior*, **36**(1), 3–15. ISSN 08943796,
[doi:10.1002/job.1960](https://doi.org/10.1002/job.1960) .

## Examples

``` r
etasqr.gen.2way(12.3, 15.6, 5.2, 7.9)
#>                                            A         B        AB
#> A treatment, B classification:      0.300000 0.5435540 0.1811847
#> A classification, B treatment:      0.484252 0.3804878 0.2047244
#> A classification, B classification: 0.300000 0.3804878 0.1268293

# Should return:
#                                           A         B        AB
# A treatment, B classification:      0.300000 0.5435540 0.1811847
# A classification, B treatment:      0.484252 0.3804878 0.2047244
# A classification, B classification: 0.300000 0.3804878 0.1268293
 
```
