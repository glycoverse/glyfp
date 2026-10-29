#' Encode Complete Glycans as Fingerprints
#'
#' @description
#' Returns one named vector per glycan using a fixed built-in motif dictionary
#' or a pretrained model. Each glycan's coordinates are independent of other
#' glycans in the input batch.
#' @param glycans A [glyrepr::glycan_structure()] vector or a character vector
#'   in any format supported by [glyparse::auto_parse()]. Formats are detected
#'   automatically, as in glymotif. For explicit format selection, use the
#'   corresponding glyparse parser first. Structures must be complete: concrete
#'   residues, specified linkages, positions, substituents and reducing-end
#'   alpha/beta state. Unknown or ambiguous fields raise an error.
#' @param method One of `"substructure_binary"` (default), `"substructure_count"`,
#'   `"sweetnet"`, `"glycangt"`, or `"gifflar"`.
#' @param ... Reserved for future options. Must be empty.
#' @param max_size Maximum residues per connected motif, an integer from 1 to 7.
#'   The default 3 selects all mono-, di- and trisaccharides, including branches.
#'   Applies only to substructure methods.
#' @param unknown For substructures, `"ignore"` (default), `"warn"` or `"error"`
#'   controls motifs absent from the fixed dictionary. Known motifs are still
#'   counted; new coordinates are never added. This does not permit uncertain
#'   input structures. Embedding methods default to and require `"error"` for
#'   tokens outside their fixed vocabularies.
#' @param na_action `"propagate"` (default) returns an all-NA vector for each
#'   missing input; `"error"` rejects missing inputs.
#' @param checkpoint For SweetNet, `"sweetnet-v2"` (default), `"sweetnet-v1"`, or
#'   the local path to either supported official artifact. Named models download
#'   from a pinned official revision on first use and are SHA-256 checked before
#'   caching and loading. GlycanGT and GIFFLAR require an explicit local file:
#'   `pretrained_tokengt_large_multitask_mask_35_final.pt` or
#'   `glycowork_gifflar_species.pt`, respectively. Obtain these yourself from
#'   [Akikitani295/GlycanGT](https://huggingface.co/Akikitani295/GlycanGT) or
#'   [DBojar/glycowork_models](https://huggingface.co/DBojar/glycowork_models).
#'   Local files are identified by SHA-256, regardless of filename, and are
#'   never replaced or downloaded. Applies only to embedding methods.
#' @param cache_dir SweetNet's named-checkpoint cache directory. No directory is
#'   created for local checkpoints or offline cache misses. Unsupported for
#'   GlycanGT and GIFFLAR.
#' @param download Whether SweetNet may download a missing named checkpoint.
#'   Set `FALSE` for offline use. GlycanGT and GIFFLAR never download weights;
#'   `FALSE` is accepted and `TRUE` is an error. Applies only to embeddings.
#'
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
#' All methods first validate complete input structures. The pretrained models
#' do not encode reducing-end alpha/beta information: only after validation,
#' embedding preprocessing omits this known state at the normal reducing carbon.
#' Thus embeddings cannot distinguish glycans differing only in this state;
#' motif fingerprints retain it. Other reducing-end positions and unsupported
#' chemistry are rejected, not guessed. Missing inputs are distinct from unknown
#' structure fields: empty input returns `list()` and missing inputs retain the
#' full coordinate set filled with NA values.
#'
#' SweetNet v1 and v2 return 128 coordinates from `bn2(lin2(...))`, before final
#' activation, dropout and classification. The default v2 follows the current
#' upstream species model; v1 and v2 coordinate spaces must not be mixed.
#' GlycanGT returns the large mask-35 model's final 768-dimensional graph token
#' and requires at least two residues. Its orthogonal random identifiers use
#' master seed 1729, a SHA-256-derived seed per canonical projected glycan, the
#' official single-graph Gaussian QR construction and width 64. This fixed
#' policy preserves caller R and CPU torch RNG state. Upstream multi-graph ORF
#' batches can yield different outputs; this adapter does not claim to reproduce
#' published downstream performance. GIFFLAR concatenates mean atom, bond and
#' residue representations (384 coordinates), using 27 fixed residue templates.
#' It rejects unsupported templates, anomeric-to-anomeric bonds and attachments
#' that would move or extend a modification. See the article for supported labels.
#'
#' Embeddings run in R with optional torch and digest; SweetNet and GlycanGT also
#' use zip. No Python, RDKit or GlyLES runtime is used. Install torch's runtime
#' explicitly: loading glyfp and calling inference never install it automatically.
#' CPU evaluation uses float32, frozen normalization, no dropout or gradients,
#' and one glycan at a time. Batch order, duplicates and batch size do not change
#' coordinates under the same model, weights, preprocessing and runtime settings.
#' Preserve these settings and package/runtime versions with saved results.
#' Recompute compared vectors when the model, dictionary or implementation changes.
#' Cross-platform comparisons should allow floating-point tolerance (the reference
#' tests bound maximum absolute coordinate error by 1e-5).
#'
#' @returns A list of named integer vectors (substructures) or double vectors
#'   (embeddings), preserving input length, order and names. Motif coordinates
#'   are named by canonical structure; embedding names are `sweetnet_001` to
#'   `sweetnet_128`, `glycangt_001` to `glycangt_768`, or `gifflar_001` to
#'   `gifflar_384`. Zero and NA entries retain all coordinates. No model or
#'   dictionary object is returned.
#' @examples
#' a <- "Gal(b1-3)Gal(b1-3)Gal(b1-"
#' b <- c("Gal(b1-4)GlcNAc(b1-", NA_character_)
#' fp_a <- glycan_fingerprint(a, "substructure_count")
#' fp_b <- glycan_fingerprint(b, "substructure_count")
#' identical(c(fp_a, fp_b), glycan_fingerprint(c(a, b), "substructure_count"))
#' fp_a[[1]][c("Gal(b1-", "Gal(b1-3)Gal(b1-")]
#' lengths(fp_b)
#' @export
glycan_fingerprint <- function(
  glycans,
  method = "substructure_binary",
  ...,
  max_size = 3L,
  unknown = "ignore",
  na_action = "propagate",
  checkpoint = "sweetnet-v2",
  cache_dir = tools::R_user_dir("glyfp", "cache"),
  download = identical(method, "sweetnet")
) {
  rlang::check_dots_empty()
  validate_options(method, max_size, unknown, na_action)
  if (method %in% c("sweetnet", "glycangt", "gifflar") && missing(unknown)) {
    unknown <- "error"
  }
  if (method %in% c("sweetnet", "glycangt", "gifflar") && !missing(max_size)) {
    cli::cli_abort(
      "{.arg max_size} applies only to substructure methods.",
      class = "glyfp_error_argument"
    )
  }
  if (
    !method %in% c("sweetnet", "glycangt", "gifflar") &&
      (!missing(checkpoint) || !missing(cache_dir) || !missing(download))
  ) {
    cli::cli_abort(
      "{.arg checkpoint}, {.arg cache_dir}, and {.arg download} apply only to embedding methods.",
      class = "glyfp_error_argument"
    )
  }
  if (method %in% c("glycangt", "gifflar")) {
    if (missing(checkpoint)) {
      cli::cli_abort(
        "{method} requires an explicit local {.arg checkpoint}; download the official weights yourself. glyfp never downloads weights for this method.",
        class = "glyfp_error_argument"
      )
    }
    if (
      !missing(cache_dir) || (!missing(download) && !identical(download, FALSE))
    ) {
      cli::cli_abort(
        "{method} uses only a user-supplied local checkpoint. {.arg cache_dir} is unsupported and {.arg download} may only be FALSE.",
        class = "glyfp_error_argument"
      )
    }
  }
  embedding <- method %in% c("sweetnet", "glycangt", "gifflar")
  if (embedding && unknown != "error") {
    cli::cli_abort(
      "Embedding methods require {.code unknown = 'error'}.",
      class = "glyfp_error_argument"
    )
  }
  glycans <- normalize_glycans(glycans, na_action)
  graphs <- as.list(glycans)
  for (i in which(!is.na(glycans))) {
    issues <- structure_issues(graphs[[i]])
    if (length(issues)) {
      cli::cli_abort(
        c(
          "Glycan at position {i} is not a supported complete structure.",
          "i" = "{paste(issues, collapse = ', ')}. Specify all residues, linkages, positions, substituents and reducing-end alpha/beta information."
        ),
        class = "glyfp_error_unsupported"
      )
    }
  }
  if (!embedding) {
    return(motif_fingerprints(glycans, method, as.integer(max_size), unknown))
  }
  glycans <- embedding_projection(glycans)
  valid <- which(!is.na(glycans))
  unique_glycans <- unique(glycans[valid])
  values <- switch(
    method,
    sweetnet = sweetnet_run(
      sweetnet_inputs(unique_glycans),
      checkpoint,
      cache_dir = cache_dir,
      download = download
    ),
    glycangt = glycangt_backend(unique_glycans, checkpoint),
    gifflar = gifflar_backend(unique_glycans, checkpoint)
  )
  features <- switch(
    method,
    sweetnet = sweetnet_features(),
    glycangt = glycangt_features(),
    gifflar = gifflar_features()
  )
  if (
    length(values) != length(unique_glycans) ||
      any(vapply(
        values,
        function(x) length(x) != length(features) || any(!is.finite(x)),
        TRUE
      ))
  ) {
    cli::cli_abort(
      "{method} returned invalid embeddings.",
      class = "glyfp_error_backend"
    )
  }
  result <- lapply(seq_along(glycans), function(i) {
    stats::setNames(rep(NA_real_, length(features)), features)
  })
  result[valid] <- lapply(match(glycans[valid], unique_glycans), function(i) {
    stats::setNames(values[[i]], features)
  })
  names(result) <- names(glycans)
  result
}

validate_options <- function(method, max_size, unknown, na_action) {
  choice <- function(x, choices, arg) {
    if (!is.character(x) || length(x) != 1L || is.na(x) || !x %in% choices) {
      cli::cli_abort(
        "{.arg {arg}} must be one of {.val {choices}}.",
        class = "glyfp_error_argument"
      )
    }
  }
  choice(
    method,
    c(
      "substructure_binary",
      "substructure_count",
      "sweetnet",
      "glycangt",
      "gifflar"
    ),
    "method"
  )
  choice(unknown, c("ignore", "warn", "error"), "unknown")
  choice(na_action, c("propagate", "error"), "na_action")
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
