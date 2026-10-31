#' Encode Residue-Centered Circular Neighborhoods
#'
#' @description
#' Computes a fixed-size hashed fingerprint from residue-level Weisfeiler-Lehman
#' (WL) neighborhoods, without fitting a dictionary or training a model.
#' @inheritParams fp_motif_count
#' @param radius Maximum neighborhood radius in glycosidic edges, an integer
#'   from 0 to 10. The default 2 includes radii 0, 1 and 2.
#' @param n_bits Number of fixed hash bins, an integer from 1 to 1048576.
#'   The default 2048 balances vector size and folding collisions; it is not
#'   optimized for a particular prediction task.
#' @param mode `"binary"` (default) records occupied bins; `"count"` sums feature
#'   occurrences, including repeated environments and collisions in each bin.
#' @details
#' Encoding version 1 starts with the concrete residue label (including the
#' stereochemistry it represents), positioned substituents, the residue's own
#' anomer/donor carbon, and a reducing-end flag. The root retains its specified
#' alpha/beta state. Each WL round hashes the previous label with a sorted
#' multiset of neighbor messages containing the neighbor label, full glycosidic
#' linkage and donor/acceptor direction. This is not an atomic stereochemistry
#' reconstruction; information absent from glyrepr labels cannot be recovered.
#'
#' One occurrence per residue per radius 0 through `radius` is recorded,
#' including repeated environments. A count fingerprint therefore sums to
#' `number_of_residues * (radius + 1)`. The iteration number distinguishes radii;
#' equal neighborhoods at different radii are not deduplicated as in some ECFP
#' implementations. Limited radius can miss differences in long-range layout.
#' This is a residue-level WL design inspired by circular molecular fingerprints,
#' not an exact implementation of atom-level Morgan/ECFP or MAP4.
#'
#' @section Fixed hashing and collisions:
#' Fields are UTF-8 strings prefixed by their byte lengths. Full SHA-256 digests
#' encode WL labels and method-specific features. The first 32 bits of each final
#' digest, interpreted big-endian, are reduced modulo `n_bits` to select a bin.
#' Colliding features add their counts or share a binary bit. There is no random
#' seed, learned vocabulary, or batch-dependent dictionary. Distinct structures
#' can share a fingerprint because of feature loss or collisions; equality is
#' not proof of isomorphism. Larger `n_bits` reduces folding collisions but does
#' not restore missing structural information. Bin names encode method, version,
#' radius, dimension and mode. Compare vectors only under identical settings and
#' compatible package/dependency versions. Counts are occurrence counts, not
#' counts of unique subgraphs. No torch runtime or model weights are used.
#' @inheritSection fp_motif_count Complete structures
#' @returns A list of named integer vectors of length `n_bits`, preserving input
#'   length, order and names. Missing inputs produce all-NA vectors and empty
#'   input produces `list()`. A single residue contributes once per radius.
#' @references
#' Rogers and Hahn (2010). Extended-connectivity fingerprints.
#' \doi{10.1021/ci100050t}.
#' Shervashidze et al. (2011). Weisfeiler-Lehman graph kernels.
#' <https://www.jmlr.org/papers/v12/shervashidze11a.html>.
#' @seealso [fp_residue_pair()], [fp_motif_count()], [fp_motif_binary()]
#' @examples
#' x <- c("Gal(b1-4)GlcNAc(b1-", "Gal(b1-3)GlcNAc(b1-")
#' fp <- fp_circular(x)
#' lengths(fp)
#' identical(fp_circular(x[1]), fp[1])
#' sum(fp_circular(x[1], mode = "count")[[1]])
#' @export
fp_circular <- function(
  glycans,
  ...,
  radius = 2L,
  n_bits = 2048L,
  mode = "binary",
  na_action = "propagate"
) {
  rlang::check_dots_empty()
  validate_residue_options(radius, n_bits, mode)
  glycans <- complete_glycans(glycans, na_action)
  hashed_residue_fingerprints(
    glycans,
    as.integer(radius),
    as.integer(n_bits),
    mode,
    "circular"
  )
}
