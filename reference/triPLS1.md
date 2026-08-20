# Tri-PLS1: three-way PLS regressed onto a y vector

Tri-PLS1: three-way PLS regressed onto a y vector

## Usage

``` r
triPLS1(X, y, numComponents, tol = 1e-10, maxIter = 100)
```

## Arguments

- X:

  Centered tensor of independent data

- y:

  Centered dependent variable

- numComponents:

  Number of components to fit

- tol:

  Relative change in loss for the model to converge (default 1e-10).

- maxIter:

  Maximum number of iterations (default 100).

## Value

Model

## Examples

``` r
set.seed(123)
X <- array(rnorm(100 * 5 * 4), dim = c(100, 5, 4))  # Random tensor (100 samples, 5 vars, 4 vars)
y <- rnorm(100)  # Random response variable
model <- triPLS1(X, y, numComponents = 2)
```
