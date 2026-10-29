
<!-- README.md is generated from README.Rmd. Please edit that file -->

# glyfp <a href="#documentation"><img src="man/figures/logo.png" align="right" height="138" alt="glyfp logo" /></a>

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-universe
version](https://glycoverse.r-universe.dev/glyfp/badges/version)](https://glycoverse.r-universe.dev/glyfp)
[![R-CMD-check](https://github.com/glycoverse/glyfp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/glycoverse/glyfp/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/glycoverse/glyfp/graph/badge.svg)](https://app.codecov.io/gh/glycoverse/glyfp)
<!-- badges: end -->

The glyfp package represents glycans as fingerprints for comparison and
modeling. It provides binary and count fingerprints of connected glycan
substructures, reusable encoders for multiple batches, and optional
pretrained SweetNet and GlycanGT embeddings. Inputs can be glyrepr
structures or text in formats supported by glyparse.

## Installation

glyfp is under local development and is not yet available from CRAN,
R-universe, or a remote repository. Install it from this package
directory:

``` r
# install.packages("pak")
pak::repo_add(glycoverse = "https://glycoverse.r-universe.dev")
pak::local_install(".")
```

**Note:** Tips and troubleshooting for the meta-package
[glycoverse](https://github.com/glycoverse/glycoverse) also apply here:
[Installation of
glycoverse](https://github.com/glycoverse/glycoverse#installation).

## Documentation

-   🚀 Get started: [Get started with glyfp](vignettes/glyfp.Rmd)
-   📚 Reference: [Fingerprint functions and
    encoders](man/glycan_fingerprint.Rd)

These links point to the local source documentation until the package
website is published. After installation, use
`help("glycan_fingerprint", package = "glyfp")` and
`browseVignettes("glyfp")`.

## Role in `glycoverse`

glyfp is being developed as part of the [Glycoverse
ecosystem](https://github.com/glycoverse/glycoverse). It builds on
glyrepr for glycan structures, glyparse for text parsing, and glymotif
for substructure extraction and counting. Its fingerprints are intended
as inputs to glycan comparison, clustering, and predictive workflows.

## Example

``` r
library(glyfp)
```

Encode a batch of glycans as substructure counts:

``` r
glycans <- c("Gal(b1-3)Gal(b1-3)Gal(b1-", "Gal(b1-")
glycan_fingerprint(glycans, method = "substructure_count", max_size = 2)
#> [[1]]
#>          Gal(b1- Gal(b1-3)Gal(b1- 
#>                3                2 
#> 
#> [[2]]
#>          Gal(b1- Gal(b1-3)Gal(b1- 
#>                1                0
```

The repeated disaccharide is counted twice, including overlapping
occurrences. Use `method = "substructure_binary"`, the default, for
presence/absence.

For multiple batches, fit once and reuse the encoder:

``` r
encoder <- fit_fingerprint(glycans, method = "substructure_binary")
transform_fingerprint(c("Gal(b1-3)Gal(b1-", NA_character_), encoder)
#> [[1]]
#>                   Gal(b1-          Gal(b1-3)Gal(b1- Gal(b1-3)Gal(b1-3)Gal(b1- 
#>                         1                         1                         0 
#> 
#> [[2]]
#>                   Gal(b1-          Gal(b1-3)Gal(b1- Gal(b1-3)Gal(b1-3)Gal(b1- 
#>                        NA                        NA                        NA
```

Each top-level substructure call fits its own dictionary; independently
fitted batches are not directly comparable, even when their vectors have
equal lengths. The [getting-started vignette](vignettes/glyfp.Rmd) shows
both combining batches in one call and reusing an encoder, along with
missing and unknown input.
