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

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

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
