# Encode Glycans as Motif Count or Binary Fingerprints

`fp_motif_count()` counts connected motifs in each glycan.
`fp_motif_binary()` records their presence or absence. By default, both
use a fixed built-in dictionary, independent of other glycans in the
input batch.

## Usage

``` r
fp_motif_count(
  glycans,
  ...,
  max_size = 3L,
  motif_dict = "fixed",
  unknown = "ignore",
  na_action = "propagate"
)

fp_motif_binary(
  glycans,
  ...,
  max_size = 3L,
  motif_dict = "fixed",
  unknown = "ignore",
  na_action = "propagate"
)
```

## Arguments

- glycans:

  A
  [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
  vector or a character vector in any format supported by
  [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html).
  Formats are detected automatically, as in glymotif. For explicit
  format selection, use the corresponding glyparse parser first.
  Structures must be complete: concrete residues, specified linkages,
  positions, substituents and reducing-end alpha/beta state. Unknown or
  ambiguous fields raise an error.

- ...:

  Reserved for future options. Must be empty.

- max_size:

  Maximum residues per connected motif, an integer from 1 to 7. The
  default 3 selects all mono-, di- and trisaccharides, including
  branches.

- motif_dict:

  `"fixed"` (default) uses the built-in glydb dictionary; `"dynamic"`
  builds a dictionary from the non-missing glycans in this call. Dynamic
  dictionaries may differ across batches, so independently computed
  vectors cannot be assumed to have aligned coordinates.

- unknown:

  `"ignore"` (default), `"warn"` or `"error"` controls motifs absent
  from the fixed dictionary. Known motifs are still counted; new
  coordinates are never added in fixed mode. In dynamic mode, all
  extracted motifs are included, so this option has no effect. Neither
  mode permits uncertain input structures.

- na_action:

  `"propagate"` (default) returns an all-NA vector for each missing
  input; `"error"` rejects missing inputs.

## Value

A list of named integer vectors, preserving input length, order and
names. Motif coordinates are named by canonical structure. Zero and NA
entries retain all coordinates. No dictionary object is returned.

## Details

The built-in dictionary `glydb-2.12.1-complete-v1` contains 24,124
canonical connected motifs of sizes 1–7, derived from 8,211 unique
complete structures (8,416 accessions) in glydb's GlyGen 2.12.1 data.
The default maximum size 3 selects 3,173 coordinates. Features are
sorted by UTF-8 bytes, without hashing or folding. The dictionary does
not grow or change with the input batch. Different `max_size` settings
or future dictionary versions define different coordinate sets. Coverage
is limited to this source snapshot; unseen motifs follow `unknown`. See
the Get started article for provenance and exclusions.

With `motif_dict = "dynamic"`, each call constructs its own dictionary
from all motifs up to `max_size` in the supplied non-missing structures,
including motifs absent from glydb. Unique features are sorted by UTF-8
bytes, so input order and duplicates do not affect the dictionary.
Different batches may produce different dictionaries, even at equal
vector lengths; combine glycans in one call to ensure aligned dynamic
coordinates. Empty input returns
[`list()`](https://rdrr.io/r/base/list.html). An all-missing dynamic
batch has no features and returns one named zero-length integer vector
per input; retain the original missingness mask because a zero-length
vector cannot hold an NA value. A missing input in a mixed batch has
all-NA entries for that batch's dynamic dictionary.

Features use
[`glymotif::extract_motif()`](https://glycoverse.github.io/glymotif/reference/extract_motif.html)
and
[`glymotif::count_motifs()`](https://glycoverse.github.io/glymotif/reference/count_motif.html)
with substructure alignment, linkages included, strict substituent
matching and strict mode. Overlapping occurrences count separately;
symmetric mappings of the same residue set count once. Binary values are
`as.integer(count > 0)`. Residue identities, linkage labels,
substituents and reducing-end alpha/beta information are retained.
Protein/lipid attachments and 3D geometry are not represented. Alditols
and floating parts/substituents are unsupported.

## Complete structures

All fingerprint functions validate complete input structures. Missing
inputs are distinct from unknown structure fields: empty input returns
[`list()`](https://rdrr.io/r/base/list.html) and missing inputs retain
the full coordinate set filled with NA values. Recompute compared
vectors when the model, dictionary or implementation changes.

## See also

[`fp_sweetnet()`](https://glycoverse.github.io/glyfp/reference/fp_sweetnet.md),
[`fp_glycangt()`](https://glycoverse.github.io/glyfp/reference/fp_glycangt.md),
[`fp_gifflar()`](https://glycoverse.github.io/glyfp/reference/fp_gifflar.md)

## Examples

``` r
a <- "Gal(b1-3)Gal(b1-3)Gal(b1-"
b <- c("Gal(b1-4)GlcNAc(b1-", NA_character_)
fp_a <- fp_motif_count(a)
fp_b <- fp_motif_count(b)
identical(c(fp_a, fp_b), fp_motif_count(c(a, b)))
#> [1] TRUE
fp_a[[1]][c("Gal(b1-", "Gal(b1-3)Gal(b1-")]
#>          Gal(b1- Gal(b1-3)Gal(b1- 
#>                3                2 
lengths(fp_b)
#> [1] 3173 3173
identical(fp_motif_binary(a), lapply(fp_a, \(x) 1L * (x > 0L)))
#> [1] TRUE
dynamic <- fp_motif_count(a, motif_dict = "dynamic")
dynamic[[1]]
#>                   Gal(b1-          Gal(b1-3)Gal(b1- Gal(b1-3)Gal(b1-3)Gal(b1- 
#>                         3                         2                         1 
```
