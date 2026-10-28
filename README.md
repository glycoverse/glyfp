
<!-- README.md is generated from README.Rmd. Please edit that file -->

# glyfp <img src="man/figures/logo.png" align="right" height="138" alt="glyfp: glycan fingerprints" />

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

**glyfp** turns glycans into structural fingerprints for comparison and
modeling. It learns a dictionary of connected residue substructures and
encodes each glycan as a binary or count vector. It uses **glymotif**
for extraction and matching, **glyparse** for text input, and
**glyrepr** for structure representation.

-   `glycan_fingerprint()` fits and transforms one batch, returning a
    simple list.
-   `fit_fingerprint()` returns a reusable encoder containing the
    dictionary, method, and parameters.
-   `transform_fingerprint()` uses that encoder for another batch.

## Installation

This is a local development version, not a CRAN or R-universe release.
From an R session in this package directory:

``` r
# install.packages("pak")
pak::repo_add(glycoverse = "https://glycoverse.r-universe.dev")
pak::local_install(".")
```

## One batch

Text formats are auto-detected exactly as in glymotif. All formats
supported by `glyparse::auto_parse()` are accepted, including
IUPAC-condensed, IUPAC-short, GlycoCT, and WURCS. A `glyrepr_structure`
vector is also accepted. Explicit format selection can be performed with
a glyparse parser before calling glyfp.

``` r
library(glyfp)
glycans <- c(short = "Galb3Galb3Galb-", single = "Gal(b1-")
fp <- glycan_fingerprint(glycans, method = "substructure_count", max_size = 2)
fp
#> $short
#>          Gal(b1- Gal(b1-3)Gal(b1- 
#>                3                2 
#> 
#> $single
#>          Gal(b1- Gal(b1-3)Gal(b1- 
#>                1                0
```

The repeated disaccharide occurs twice, with one shared residue.
Symmetric mappings of the same residue set count once.
`substructure_binary` is the default and applies `count > 0` to the same
counts. `max_size = 3` is the default, following glymotif, and includes
all one-, two-, and three-residue fragments, including branched
trisaccharides. Values from 1 through 7 are accepted; large, branched
structures become expensive as the maximum size grows.

## Comparable fingerprints across batches

**Every `glycan_fingerprint()` call fits its own batch. Results from
independent calls cannot be compared directly, even if the vectors have
equal lengths.** Their positions can refer to different motifs.

For a joint comparison, combine all glycans in one call:

``` r
a <- c("Gal(b1-3)GalNAc(a1-", "GalNAc(a1-")
b <- "Gal(b1-4)GlcNAc(b1-"
all_fp <- glycan_fingerprint(c(a, b))
lengths(all_fp)
#> [1] 5 5 5
```

For repeated or predictive use, fit once and transform each batch with
the same encoder. Fit on training data only when evaluating a predictive
model.

``` r
encoder <- fit_fingerprint(a, method = "substructure_count")
encoder
#> <glyfp_encoder>
#> Method: substructure_count
#> Features: 3
#> Maximum residues: 3
encoder$features
#> [1] "Gal(b1-"             "Gal(b1-3)GalNAc(a1-" "GalNAc(a1-"
fp_a <- transform_fingerprint(a, encoder)
fp_b <- transform_fingerprint(b, encoder)
identical(names(fp_a[[1]]), names(fp_b[[1]]))
#> [1] TRUE
```

Save and reload the encoder with `saveRDS()` / `readRDS()`. It stores
ordinary R data and dependency versions. Transform rejects a
dependency-version mismatch; restore the recorded environment or fit a
new encoder and recompute every batch.

``` r
path <- tempfile(fileext = ".rds")
saveRDS(encoder, path)
identical(transform_fingerprint(b, readRDS(path)), fp_b)
#> [1] TRUE
unlink(path)
```

## What each feature means

Feature names are canonical IUPAC motif strings, sorted by UTF-8 bytes,
with no hashing or folding. Counts use glymotif’s substructure
alignment, with linkages included, strict substituent matching, and
strict matching mode. This preserves branch topology, residue
identities, linkage positions, and anomer constraints. A motif
half-linkage constrains an internal residue’s anomer as well as a
reducing-end residue’s anomer.

Generic motif residues match compatible concrete residues; unknown motif
linkage positions act as wildcards. Thus feature columns can overlap in
meaning. The reverse direction is not assumed: a generic target is not
evidence of a specific residue. These are glymotif matching semantics,
not atom-level or 3D descriptors.

## Missing, novel, and unsupported input

-   `NA` does not contribute features while fitting. Transformation
    preserves its position with an all-NA integer vector. Use
    `na_action = "error"` to reject it.
-   Unknown features do not grow the dictionary. `unknown = "ignore"` is
    the default; `"warn"` or `"error"` detects new extracted canonical
    motifs. A new specific motif can still match a broader motif already
    in the dictionary.
-   Empty input returns `list()`. Empty/all-NA training data produces a
    zero-feature encoder. Every output vector then has length zero, so
    retain `is.na(glycans)` if you need to distinguish missing entries
    in this special case.
-   Invalid text is an error with the original glyparse error attached.
    Names, duplicates, and input order are preserved; duplicate
    structures share the same computed values.
-   Unresolved floating parts, floating substituents, and alditols are
    explicitly rejected in this first version. Embedding methods are not
    yet part of the R API; an actual pretrained SweetNet CPU feasibility
    probe is included under `dev/`, with its environment and checkpoint
    recorded in the method review.

See [the method review](dev/methods.md) for literature, scope, and the R
torch versus Python implementation assessment. The [fingerprint
vignette](vignettes/fingerprints.Rmd) contains further executable
examples.

## Development

The package layout, generated README, pkgdown index, and R-CMD-check
workflow follow the development package
[glyenzy](https://github.com/glycoverse/glyenzy). CI checks Linux,
macOS, and Windows plus R devel/oldrel and saves local coverage as a
workflow artifact. It has no publishing or release step. Remote badges
will be added after a repository exists and its runs can be verified.

``` r
devtools::document()
devtools::test()
devtools::build_readme()
devtools::check()
```
