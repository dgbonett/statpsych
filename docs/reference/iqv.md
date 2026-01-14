# Indices of qualitative variation

Computes the Shannon, Berger, and Simpson indices of qualitative
variation. Will be replaced soon with ci.diversity.

## Usage

``` r
iqv(f)
```

## Arguments

- f:

  vector of multinomial frequency counts

## Value

Returns estimates of the Shannon, Berger, and Simpson indices

## Examples

``` r
f <- c(10, 46, 15, 3)
iqv(f)
#>    Simpson    Berger   Shannon
#>  0.7367908 0.5045045 0.7353931

# Should return:
#   Simpson    Berger   Shannon
# 0.7367908 0.5045045       0.7
 
```
