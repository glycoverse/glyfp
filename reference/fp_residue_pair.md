# Encode Residue Environments and Graph Distances

Computes a fixed-size hashed fingerprint of unordered residue pairs
using their endpoint environments and shortest glycosidic-edge distance.

## Usage

``` r
fp_residue_pair(
  glycans,
  ...,
  radius = 1L,
  n_bits = 2048L,
  mode = "binary",
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

  Radius of each endpoint's WL environment, an integer from 0 to 10. The
  default 1 captures immediate neighbors and linkage labels. Only the
  final radius is used for endpoints, not all intermediate radii.

- n_bits:

  Number of fixed hash bins, an integer from 1 to 1048576. The default
  2048 balances vector size and folding collisions; it is not optimized
  for a particular prediction task.

- mode:

  `"binary"` (default) records occupied bins; `"count"` sums feature
  occurrences, including repeated environments and collisions in each
  bin.

- na_action:

  `"propagate"` (default) returns an all-NA vector for each missing
  input; `"error"` rejects missing inputs.

## Value

A list of named integer vectors of length `n_bits`, preserving input
length, order and names. Missing inputs yield all-NA vectors, a single
residue yields all zeros, and empty input yields
[`list()`](https://rdrr.io/r/base/list.html).

## Details

Every distinct unordered residue pair contributes one feature consisting
of the lexically sorted endpoint labels and their undirected
shortest-path length. Endpoints use the same residue labels, root
information and directed linkage messages as
[`fp_circular_binary()`](https://glycoverse.github.io/glyfp/reference/fp_circular_binary.md).
All distances are included; there is no cutoff. Distance is a number of
glycosidic edges, not a 3D or spatial distance.

A count vector sums to `n * (n - 1) / 2` for `n` residues. There are no
self-pairs: a single-residue glycan yields an all-zero vector and cannot
be distinguished from another single residue with this method. Use
[`fp_circular_binary()`](https://glycoverse.github.io/glyfp/reference/fp_circular_binary.md)
or motif features when that distinction matters.

Endpoint environments can preserve local linkage and anomer information,
but the feature does not record the sequence of residues or linkages
along the intervening path beyond those environments. With `radius = 0`,
even linkage acceptor positions are omitted. Pair histograms are not
complete graph invariants. Large glycans require quadratic pair
enumeration and an all-pairs distance matrix. This residue-level design
is inspired by atom-pair features; it is not an implementation of
atom-level atom pairs or MAP4.

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

Carhart, Smith and Venkataraghavan (1985). Atom pairs as molecular
features in structure-activity studies: definition and applications.
[doi:10.1021/ci00046a002](https://doi.org/10.1021/ci00046a002) .

## See also

[`fp_circular_binary()`](https://glycoverse.github.io/glyfp/reference/fp_circular_binary.md),
[`fp_motif_count()`](https://glycoverse.github.io/glyfp/reference/fp_motif_count.md),
[`fp_motif_binary()`](https://glycoverse.github.io/glyfp/reference/fp_motif_count.md)

## Examples

``` r
x <- c("Gal(b1-4)GlcNAc(b1-", "Gal(b1-4)GlcNAc(b1-4)GlcNAc(b1-")
fp <- fp_residue_pair(x, mode = "count")
vapply(fp, sum, 0)
#> [1] 1 3
identical(fp_residue_pair(x[1], mode = "count"), fp[1])
#> [1] TRUE
```
