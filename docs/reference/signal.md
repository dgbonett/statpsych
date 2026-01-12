# Parameter estimates for a signal detection study

Computes the hit rate, false alarm rate, d-prime, threshold, and bias
for one participant (observer) in a Yes/No signal detection study. An
equal-variance Gaussian model is assumed. The parameter estimates are
computed after adding .5 to the number of "Yes" responses in each
condition (the signal and noise conditions) and adding 1 to the number
of signal trials and to the number of noise trails. In memory
recognition studies, the observer is first presented with set of words
or images to study, and is later presented with another set of words or
images where some items are from the first list (old items) and some
items are new items.

For more details, see Section 3.8 of Bonett (2021, Volume 3)

## Usage

``` r
signal(f1, f2, n1, n2)
```

## Arguments

- f1:

  number of "Yes" responses in the stimulus (old item) trials

- f2:

  number of "Yes" responses in the noise (new item) trials

- n1:

  number of stimulus (or old item) trials

- n2:

  number of noise (or new item) trials

## Value

Returns a 1-row matrix. The columns are:

- HR - estimate of hit rate

- FAR - estimate of false alarm rate

- d-prime - estimate of d-prime

- Threshold - estimate of threshold (criterion)

- Bias - estimate of threshold minus d-prime/2

## References

Wickens TD (2002). *Elementary Signal Detection Theory*. Oxford.

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
signal(82, 46, 100, 100)
#>         HR      FAR  d-prime  Threshold       Bias
#>  0.8168317 0.460396 1.002793 0.09943603 -0.4019603

# Should return:
#         HR      FAR  d-prime  Threshold       Bias
#  0.8168317 0.460396 1.002793 0.09943603 -0.4019603

```
