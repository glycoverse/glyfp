# Encode Glycans with GlycanGT

Returns 768-dimensional graph Transformer embeddings from the official
GlycanGT large mask-35 model using a manually supplied local checkpoint.

## Usage

``` r
fp_glycangt(glycans, checkpoint, ..., na_action = "propagate")
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

  Local path to the official
  `pretrained_tokengt_large_multitask_mask_35_final.pt` artifact. Obtain
  it from
  [Akikitani295/GlycanGT](https://huggingface.co/Akikitani295/GlycanGT).
  The file is identified by SHA-256, regardless of filename, and is
  never replaced or downloaded. A local checkpoint is required on every
  call.

- ...:

  Reserved for future options. Must be empty.

- na_action:

  `"propagate"` (default) returns an all-NA vector for each missing
  input; `"error"` rejects missing inputs.

## Value

A list of named double vectors with 768 coordinates (`glycangt_001` to
`glycangt_768`), preserving input length, order and names. Missing
inputs return all-NA vectors with the full coordinate set. No model
object is returned.

## Details

GlycanGT returns the model's final graph token and requires at least two
residues. Its orthogonal random identifiers use master seed 1729, a
SHA-256-derived seed per canonical projected glycan, the official
single-graph Gaussian QR construction and width 64. This fixed policy
preserves caller R and CPU torch RNG state. Upstream multi-graph ORF
batches can yield different outputs; this adapter does not claim to
reproduce published downstream performance.

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
[`fp_gifflar()`](https://glycoverse.github.io/glyfp/reference/fp_gifflar.md)

## Examples

``` r
checkpoint <- Sys.getenv("GLYFP_TEST_GLYCANGT_CHECKPOINT")
if (file.exists(checkpoint) && requireNamespace("torch", quietly = TRUE) &&
    torch::torch_is_installed()) {
  fp <- fp_glycangt("Gal(b1-4)GlcNAc(b1-", checkpoint = checkpoint)
  lengths(fp)
}
```
