# NPLStoolbox

## Overview

The `NPLStoolbox` allows researchers to use the N-way Partial Least
Squares method for their multi-way data.

- [`ncrossreg()`](http://grvanderploeg.com/NPLStoolbox/reference/ncrossreg.md)
  allows the user to identify the appropriate number of NPLS components
  for their data.
- [`triPLS1()`](http://grvanderploeg.com/NPLStoolbox/reference/triPLS1.md)
  allows the user to create an NPLS model.
- [`npred()`](http://grvanderploeg.com/NPLStoolbox/reference/npred.md)
  allows the user to predict y for new data.

This package also comes with two example datasets:

- `Cornejo2025`: a clinical observational cohort study of 39 transgender
  persons starting gender-affirming hormone therapy, containing
  longitudinally measured tongue microbiome, salivary microbiome,
  salivary cytokine, salivary biochemistry, and circulatory hormone
  levels (doi TBD).
- `Jakobsen2025`: an observational cohort of 169 mother-infant dyads
  investigating the effect of maternal obesity on human milk and the
  infant gut microbiome <https://doi.org/10.21203/rs.3.rs-6244750/v1>.

## Documentation

A basic introduction to the package using the example dataset is given
in
[`vignette("Introduction")`](http://grvanderploeg.com/NPLStoolbox/articles/Introduction.md).

This vignette and all function documentation can be found
[here](https://grvanderploeg.com/NPLStoolbox/).

## Installation

The `NPLStoolbox` package can be installed from CRAN using:

``` r

install.packages("NPLStoolbox")
```

## Development version

You can install the development version of NPLStoolbox from
[GitHub](https://github.com/) with:

``` r

# install.packages("pak")
pak::pak("GRvanderPloeg/NPLStoolbox")
```

## Usage

``` r

library(parafac4microbiome)
library(NPLStoolbox)
set.seed(123)

# Process one of the data cubes from Cornejo2025
processedTongue = processDataCube(Cornejo2025$Tongue_microbiome, sparsityThreshold=0.5, considerGroups=TRUE, groupVariable="GenderID", centerMode=1, scaleMode=2)

# Prepare Y: binarized gender identity
Y = as.numeric(as.factor(Cornejo2025$Tongue_microbiome$mode1$GenderID))
Ycnt = Y - mean(Y)

# Make a one-component NPLS model
model = triPLS1(processedTongue$data, Ycnt, 1)
```

## Getting help

If you encounter an unexpected error or a clear bug, please file an
issue with a minimal reproducible example here on
[Github](https://github.com/GRvanderPloeg/NPLStoolbox/issues). For
questions or other types of feedback, feel free to send an email.
