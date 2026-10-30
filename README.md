
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
substructures and optional pretrained SweetNet, GlycanGT and GIFFLAR
embeddings through `fp_motif_binary()`, `fp_motif_count()`,
`fp_sweetnet()`, `fp_glycangt()` and `fp_gifflar()`. Each function
exposes only options relevant to its method. Inputs can be glyrepr
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
-   📚 Reference: [Motif fingerprints](man/fp_motif_count.Rd),
    [SweetNet](man/fp_sweetnet.Rd), [GlycanGT](man/fp_glycangt.Rd),
    [GIFFLAR](man/fp_gifflar.Rd)

These links point to the local source documentation until the package
website is published. After installation, use
`help("fp_motif_count", package = "glyfp")` and
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

Encode complete glycans using the built-in glydb motif dictionary:

``` r
glycans <- c("Gal(b1-3)Gal(b1-3)Gal(b1-", "Gal(b1-")
fp <- fp_motif_count(glycans)
lengths(fp)
#> [1] 3173 3173
fp[[1]][c("Gal(b1-", "Gal(b1-3)Gal(b1-")]
#>          Gal(b1- Gal(b1-3)Gal(b1- 
#>                3                2
```

The repeated disaccharide is counted twice, including overlapping
occurrences. Use `fp_motif_binary(glycans)` for presence/absence. The
default `motif_dict = "fixed"` keeps coordinates aligned across separate
calls. Set `motif_dict = "dynamic"` in either motif function to build a
dictionary from the current batch, including motifs outside glydb.
Dynamic dictionaries may differ between batches; combine inputs in one
call for aligned coordinates. Inputs must specify complete structures,
including the reducing-end anomer. The [getting-started
vignette](vignettes/glyfp.Rmd) explains all methods, pretrained
checkpoints and supported chemistry.
