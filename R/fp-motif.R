#' Encode Glycans as Motif Count or Binary Fingerprints
#'
#' @description
#' `fp_motif_count()` counts connected motifs in each glycan using a fixed
#' built-in dictionary. `fp_motif_binary()` records their presence or absence.
#' Each glycan's coordinates are independent of other glycans in the input batch.
#' @param glycans A [glyrepr::glycan_structure()] vector or a character vector
#'   in any format supported by [glyparse::auto_parse()]. Formats are detected
#'   automatically, as in glymotif. For explicit format selection, use the
#'   corresponding glyparse parser first. Structures must be complete: concrete
#'   residues, specified linkages, positions, substituents and reducing-end
#'   alpha/beta state. Unknown or ambiguous fields raise an error.
#' @param ... Reserved for future options. Must be empty.
#' @param max_size Maximum residues per connected motif, an integer from 1 to 7.
#'   The default 3 selects all mono-, di- and trisaccharides, including branches.
#' @param unknown `"ignore"` (default), `"warn"` or `"error"` controls motifs
#'   absent from the fixed dictionary. Known motifs are still counted; new
#'   coordinates are never added. This does not permit uncertain input structures.
#' @param na_action `"propagate"` (default) returns an all-NA vector for each
#'   missing input; `"error"` rejects missing inputs.
#' @details
#' The built-in dictionary `glydb-2.12.1-complete-v1` contains 24,124 canonical
#' connected motifs of sizes 1--7, derived from 8,211 unique complete structures
#' (8,416 accessions) in glydb's GlyGen 2.12.1 data. The default maximum size 3
#' selects 3,173 coordinates. Features are sorted by UTF-8 bytes, without hashing
#' or folding. The dictionary does not grow or change with the input batch.
#' Different `max_size` settings or future dictionary versions define different
#' coordinate sets. Coverage is limited to this source snapshot; unseen motifs
#' follow `unknown`. See the Get started article for provenance and exclusions.
#'
#' Features use [glymotif::extract_motif()] and [glymotif::count_motifs()] with
#' substructure alignment, linkages included, strict substituent matching and
#' strict mode. Overlapping occurrences count separately; symmetric mappings of
#' the same residue set count once. Binary values are `as.integer(count > 0)`.
#' Residue identities, linkage labels, substituents and reducing-end alpha/beta
#' information are retained. Protein/lipid attachments and 3D geometry are not
#' represented. Alditols and floating parts/substituents are unsupported.
#'
#' @section Complete structures:
#' All fingerprint functions validate complete input structures. Missing inputs
#' are distinct from unknown structure fields: empty input returns `list()` and
#' missing inputs retain the full coordinate set filled with NA values.
#' Recompute compared vectors when the model, dictionary or implementation changes.
#'
#' @returns A list of named integer vectors, preserving input length, order and
#'   names. Motif coordinates are named by canonical structure. Zero and NA
#'   entries retain all coordinates. No dictionary object is returned.
#' @seealso [fp_sweetnet()], [fp_glycangt()], [fp_gifflar()]
#' @examples
#' a <- "Gal(b1-3)Gal(b1-3)Gal(b1-"
#' b <- c("Gal(b1-4)GlcNAc(b1-", NA_character_)
#' fp_a <- fp_motif_count(a)
#' fp_b <- fp_motif_count(b)
#' identical(c(fp_a, fp_b), fp_motif_count(c(a, b)))
#' fp_a[[1]][c("Gal(b1-", "Gal(b1-3)Gal(b1-")]
#' lengths(fp_b)
#' identical(fp_motif_binary(a), lapply(fp_a, \(x) 1L * (x > 0L)))
#' @export
fp_motif_count <- function(
  glycans,
  ...,
  max_size = 3L,
  unknown = "ignore",
  na_action = "propagate"
) {
  rlang::check_dots_empty()
  validate_motif_options(max_size, unknown)
  glycans <- complete_glycans(glycans, na_action)
  motif_fingerprints(glycans, as.integer(max_size), unknown, binary = FALSE)
}

#' @rdname fp_motif_count
#' @export
fp_motif_binary <- function(
  glycans,
  ...,
  max_size = 3L,
  unknown = "ignore",
  na_action = "propagate"
) {
  rlang::check_dots_empty()
  validate_motif_options(max_size, unknown)
  glycans <- complete_glycans(glycans, na_action)
  motif_fingerprints(glycans, as.integer(max_size), unknown, binary = TRUE)
}

validate_motif_options <- function(max_size, unknown) {
  validate_choice(unknown, c("ignore", "warn", "error"), "unknown")
  if (
    !is.numeric(max_size) ||
      length(max_size) != 1L ||
      is.na(max_size) ||
      !is.finite(max_size) ||
      max_size < 1 ||
      max_size > 7 ||
      max_size != as.integer(max_size)
  ) {
    cli::cli_abort(
      "{.arg max_size} must be an integer from 1 to 7.",
      class = "glyfp_error_argument"
    )
  }
}
