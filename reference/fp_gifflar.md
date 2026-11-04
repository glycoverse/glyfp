# Encode Glycans with GIFFLAR

Returns 384-dimensional molecular graph embeddings from the official
GIFFLAR species model with verified checkpoint downloads and caching.

## Usage

``` r
fp_gifflar(
  glycans,
  checkpoint = "gifflar-species",
  ...,
  cache_dir = tools::R_user_dir("glyfp", "cache"),
  download = TRUE,
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

- checkpoint:

  `"gifflar-species"` (default), or the local path to the official
  `glycowork_gifflar_species.pt` artifact. The named model downloads
  from a pinned official revision on first use and is SHA-256 checked
  before caching and loading. Local files are identified by SHA-256,
  regardless of filename, and are never replaced or downloaded.

- ...:

  Reserved for future options. Must be empty.

- cache_dir:

  Named-checkpoint cache directory. No directory is created for local
  checkpoints or offline cache misses.

- download:

  Whether to download a missing named checkpoint. Set `FALSE` for
  offline use.

- na_action:

  `"propagate"` (default) returns an all-NA vector for each missing
  input; `"error"` rejects missing inputs.

## Value

A list of named double vectors with 384 coordinates (`gifflar_001` to
`gifflar_384`), preserving input length, order and names. Missing inputs
return all-NA vectors with the full coordinate set. No model object is
returned.

## Details

GIFFLAR concatenates mean atom, bond and residue representations, using
27 fixed residue templates. It rejects unsupported templates,
anomeric-to-anomeric bonds and attachments that would move or extend a
modification. See the Get started article for supported labels.

## Embedding contract

Pretrained models do not encode reducing-end alpha/beta information:
only after complete-structure validation, preprocessing omits this known
state at the normal reducing carbon. Thus embeddings cannot distinguish
glycans differing only in this state; motif fingerprints retain it.
Other reducing-end positions and unsupported chemistry are rejected.
Tokens outside the fixed model vocabulary always raise an error.

Embeddings run in R with optional torch and digest; SweetNet and
GlycanGT also use zip. No Python, RDKit or GlyLES runtime is used.
Install torch's runtime explicitly: loading glyfp and calling inference
never install it automatically. CPU evaluation uses float32, frozen
normalization, no dropout or gradients, and one glycan at a time. Batch
order, duplicates and batch size do not change coordinates under the
same model, weights, preprocessing and runtime settings. Preserve these
settings and package/runtime versions with saved results. Cross-platform
comparisons should allow floating-point tolerance (the reference tests
bound maximum absolute coordinate error by 1e-5).

## Complete structures

All fingerprint functions validate complete input structures. Missing
inputs are distinct from unknown structure fields: empty input returns
[`list()`](https://rdrr.io/r/base/list.html) and missing inputs retain
the full coordinate set filled with NA values. Recompute compared
vectors when the model, dictionary or implementation changes.

## See also

[`fp_motif_count()`](https://glycoverse.github.io/glyfp/reference/fp_motif_count.md),
[`fp_sweetnet()`](https://glycoverse.github.io/glyfp/reference/fp_sweetnet.md),
[`fp_glycangt()`](https://glycoverse.github.io/glyfp/reference/fp_glycangt.md)

## Examples

``` r
checkpoint <- Sys.getenv("GLYFP_TEST_GIFFLAR_CHECKPOINT")
if (file.exists(checkpoint) && requireNamespace("torch", quietly = TRUE) &&
    torch::torch_is_installed()) {
  fp <- fp_gifflar("Gal6S(b1-4)GlcNAc(b1-", checkpoint = checkpoint)
  lengths(fp)
}
```
