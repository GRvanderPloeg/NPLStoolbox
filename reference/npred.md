# Predict Y for new data by projecting the data onto the latent space defined by an NPLS model.

Predict Y for new data by projecting the data onto the latent space
defined by an NPLS model.

## Usage

``` r
npred(model, newX)
```

## Arguments

- model:

  NPLS model

- newX:

  New data organized in a matrix of (Inew x J x K) with Inew new
  subjects

## Value

Ypred: vector of the predicted value(s) of Y for the new data

## Examples

``` r
Y = as.numeric(as.factor(Cornejo2025$Tongue$mode1$GenderID))
Ycnt = Y - mean(Y)
model = triPLS1(Cornejo2025$Tongue$data, Ycnt, numComponents=1)
npred(model, Cornejo2025$Tongue$data[1,,])
#>             [,1]
#> [1,] -0.00637789
```
