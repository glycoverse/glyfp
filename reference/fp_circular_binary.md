# Encode Residue-Centered Circular Neighborhoods

Computes a fixed-size hashed fingerprint from residue-level
Weisfeiler-Lehman (WL) neighborhoods, without fitting a dictionary or
training a model. `fp_circular_binary()` records occupied bins;
`fp_circular_count()` sums feature occurrences, including repeated
environments and collisions.

## Usage

``` r
fp_circular_binary(
  glycans,
  ...,
  radius = 2L,
  n_bits = 2048L,
  na_action = "propagate"
)

fp_circular_count(
  glycans,
  ...,
  radius = 2L,
  n_bits = 2048L,
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

- radius:

  Maximum neighborhood radius in glycosidic edges, an integer from 0
  to 10. The default 2 includes radii 0, 1 and 2.

- n_bits:

  Number of fixed hash bins, an integer from 1 to 1048576. The default
  2048 balances vector size and folding collisions; it is not optimized
  for a particular prediction task.

- na_action:

  `"propagate"` (default) returns an all-NA vector for each missing
  input; `"error"` rejects missing inputs.

## Value

A list of named integer vectors of length `n_bits`, preserving input
length, order and names. Missing inputs produce all-NA vectors and empty
input produces [`list()`](https://rdrr.io/r/base/list.html). A single
residue contributes once per radius.

## Details

Encoding version 1 starts with the concrete residue label (including the
stereochemistry it represents), positioned substituents, the residue's
own anomer/donor carbon, and a reducing-end flag. The root retains its
specified alpha/beta state. Each WL round hashes the previous label with
a sorted multiset of neighbor messages containing the neighbor label,
full glycosidic linkage and donor/acceptor direction. This is not an
atomic stereochemistry reconstruction; information absent from glyrepr
labels cannot be recovered.

One occurrence per residue per radius 0 through `radius` is recorded,
including repeated environments. A count fingerprint therefore sums to
`number_of_residues * (radius + 1)`. The iteration number distinguishes
radii; equal neighborhoods at different radii are not deduplicated as in
some ECFP implementations. Limited radius can miss differences in
long-range layout. This is a residue-level WL design inspired by
circular molecular fingerprints, not an exact implementation of
atom-level Morgan/ECFP or MAP4.

## Fixed hashing and collisions

Fields are UTF-8 strings prefixed by their byte lengths. Full SHA-256
digests encode WL labels and method-specific features. The first 32 bits
of each final digest, interpreted big-endian, are reduced modulo
`n_bits` to select a bin. Colliding features add their counts or share a
binary bit. There is no random seed, learned vocabulary, or
batch-dependent dictionary. Distinct structures can share a fingerprint
because of feature loss or collisions; equality is not proof of
isomorphism. Larger `n_bits` reduces folding collisions but does not
restore missing structural information. Bin names encode method,
version, radius, dimension and mode. Compare vectors only under
identical settings and compatible package/dependency versions. Counts
are occurrence counts, not counts of unique subgraphs. No torch runtime
or model weights are used.

## Complete structures

All fingerprint functions validate complete input structures. Missing
inputs are distinct from unknown structure fields: empty input returns
[`list()`](https://rdrr.io/r/base/list.html) and missing inputs retain
the full coordinate set filled with NA values. Recompute compared
vectors when the model, dictionary or implementation changes.

## References

Rogers and Hahn (2010). Extended-connectivity fingerprints.
[doi:10.1021/ci100050t](https://doi.org/10.1021/ci100050t) .
Shervashidze et al. (2011). Weisfeiler-Lehman graph kernels.
<https://www.jmlr.org/papers/v12/shervashidze11a.html>.

## See also

[`fp_residue_pair()`](https://glycoverse.github.io/glyfp/reference/fp_residue_pair.md),
[`fp_motif_count()`](https://glycoverse.github.io/glyfp/reference/fp_motif_count.md),
[`fp_motif_binary()`](https://glycoverse.github.io/glyfp/reference/fp_motif_count.md)

## Examples

``` r
x <- c("Gal(b1-4)GlcNAc(b1-", "Gal(b1-3)GlcNAc(b1-")
fp <- fp_circular_binary(x)
lengths(fp)
#> [1] 2048 2048
identical(fp_circular_binary(x[1]), fp[1])
#> [1] TRUE
sum(fp_circular_count(x[1])[[1]])
#> [1] 6
```
