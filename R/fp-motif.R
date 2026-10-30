#' Encode Glycans as Motif Count or Binary Fingerprints
#'
#' @description
#' `fp_motif_count()` counts connected motifs in each glycan.
#' `fp_motif_binary()` records their presence or absence. By default, both use
#' a fixed built-in dictionary, independent of other glycans in the input batch.
#' @param glycans A [glyrepr::glycan_structure()] vector or a character vector
#'   in any format supported by [glyparse::auto_parse()]. Formats are detected
#'   automatically, as in glymotif. For explicit format selection, use the
#'   corresponding glyparse parser first. Structures must be complete: concrete
#'   residues, specified linkages, positions, substituents and reducing-end
#'   alpha/beta state. Unknown or ambiguous fields raise an error.
#' @param ... Reserved for future options. Must be empty.
#' @param max_size Maximum residues per connected motif, an integer from 1 to 7.
#'   The default 3 selects all mono-, di- and trisaccharides, including branches.
#' @param motif_dict `"fixed"` (default) uses the built-in glydb dictionary;
#'   `"dynamic"` builds a dictionary from the non-missing glycans in this call.
#'   Dynamic dictionaries may differ across batches, so independently computed
#'   vectors cannot be assumed to have aligned coordinates.
#' @param unknown `"ignore"` (default), `"warn"` or `"error"` controls motifs
#'   absent from the fixed dictionary. Known motifs are still counted; new
#'   coordinates are never added in fixed mode. In dynamic mode, all extracted
#'   motifs are included, so this option has no effect. Neither mode permits
#'   uncertain input structures.
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
#' With `motif_dict = "dynamic"`, each call constructs its own dictionary from
#' all motifs up to `max_size` in the supplied non-missing structures, including
#' motifs absent from glydb. Unique features are sorted by UTF-8 bytes, so input
#' order and duplicates do not affect the dictionary. Different batches may
#' produce different dictionaries, even at equal vector lengths; combine glycans
#' in one call to ensure aligned dynamic coordinates. Empty input returns
#' `list()`. An all-missing dynamic batch has no features and returns one named
#' zero-length integer vector per input; retain the original missingness mask
#' because a zero-length vector cannot hold an NA value. A missing input in a
#' mixed batch has all-NA entries for that batch's dynamic dictionary.
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
#' dynamic <- fp_motif_count(a, motif_dict = "dynamic")
#' dynamic[[1]]
#' @export
fp_motif_count <- function(
  glycans,
  ...,
  max_size = 3L,
  motif_dict = "fixed",
  unknown = "ignore",
  na_action = "propagate"
) {
  rlang::check_dots_empty()
  validate_motif_options(max_size, unknown, motif_dict)
  glycans <- complete_glycans(glycans, na_action)
  motif_fingerprints(
    glycans,
    as.integer(max_size),
    unknown,
    binary = FALSE,
    motif_dict = motif_dict
  )
}

#' @rdname fp_motif_count
#' @export
fp_motif_binary <- function(
  glycans,
  ...,
  max_size = 3L,
  motif_dict = "fixed",
  unknown = "ignore",
  na_action = "propagate"
) {
  rlang::check_dots_empty()
  validate_motif_options(max_size, unknown, motif_dict)
  glycans <- complete_glycans(glycans, na_action)
  motif_fingerprints(
    glycans,
    as.integer(max_size),
    unknown,
    binary = TRUE,
    motif_dict = motif_dict
  )
}

validate_motif_options <- function(max_size, unknown, motif_dict) {
  validate_choice(motif_dict, c("fixed", "dynamic"), "motif_dict")
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
