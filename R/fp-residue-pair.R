#' Encode Residue Environments and Graph Distances
#'
#' @description
#' Computes a fixed-size hashed fingerprint of unordered residue pairs using
#' their endpoint environments and shortest glycosidic-edge distance.
#' @inheritParams fp_circular_binary
#' @param mode `"binary"` (default) records occupied bins; `"count"` sums feature
#'   occurrences, including repeated environments and collisions in each bin.
#' @param radius Radius of each endpoint's WL environment, an integer from 0 to
#'   10. The default 1 captures immediate neighbors and linkage labels. Only the
#'   final radius is used for endpoints, not all intermediate radii.
#' @details
#' Every distinct unordered residue pair contributes one feature consisting of
#' the lexically sorted endpoint labels and their undirected shortest-path length.
#' Endpoints use the same residue labels, root information and directed linkage
#' messages as [fp_circular_binary()]. All distances are included; there is no cutoff.
#' Distance is a number of glycosidic edges, not a 3D or spatial distance.
#'
#' A count vector sums to `n * (n - 1) / 2` for `n` residues. There are no
#' self-pairs: a single-residue glycan yields an all-zero vector and cannot be
#' distinguished from another single residue with this method. Use
#' [fp_circular_binary()] or motif features when that distinction matters.
#'
#' Endpoint environments can preserve local linkage and anomer information, but
#' the feature does not record the sequence of residues or linkages along the
#' intervening path beyond those environments. With `radius = 0`, even linkage
#' acceptor positions are omitted. Pair histograms are not complete graph
#' invariants. Large glycans require quadratic pair enumeration and an all-pairs
#' distance matrix. This residue-level design is inspired by atom-pair features;
#' it is not an implementation of atom-level atom pairs or MAP4.
#' @inheritSection fp_circular_binary Fixed hashing and collisions
#' @inheritSection fp_motif_count Complete structures
#' @returns A list of named integer vectors of length `n_bits`, preserving input
#'   length, order and names. Missing inputs yield all-NA vectors, a single
#'   residue yields all zeros, and empty input yields `list()`.
#' @references
#' Carhart, Smith and Venkataraghavan (1985). Atom pairs as molecular features in
#' structure-activity studies: definition and applications.
#' \doi{10.1021/ci00046a002}.
#' @seealso [fp_circular_binary()], [fp_motif_count()], [fp_motif_binary()]
#' @examples
#' x <- c("Gal(b1-4)GlcNAc(b1-", "Gal(b1-4)GlcNAc(b1-4)GlcNAc(b1-")
#' fp <- fp_residue_pair(x, mode = "count")
#' vapply(fp, sum, 0)
#' identical(fp_residue_pair(x[1], mode = "count"), fp[1])
#' @export
fp_residue_pair <- function(
  glycans,
  ...,
  radius = 1L,
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
    "residue_pair"
  )
}
