#' Encode Glycans as Structural Fingerprints
#'
#' @description
#' `glycan_fingerprint()` fits a dictionary on the supplied batch, then transforms
#' that batch, returning only a list of vectors. `fit_fingerprint()` creates a
#' reusable encoder; `transform_fingerprint()` applies its fixed dictionary,
#' method, and parameters to another batch.
#'
#' @param glycans A [glyrepr::glycan_structure()] vector (class
#'   `glyrepr_structure`) or a character vector in any format supported by
#'   [glyparse::auto_parse()]. Like glymotif, text formats are auto-detected;
#'   there is no separate format argument. For explicit format selection, call
#'   the corresponding glyparse parser first. Invalid text raises a parse error.
#' @param method Either `"substructure_binary"` (default, presence/absence) or
#'   `"substructure_count"` (number of occurrences).
#' @param ... For `glycan_fingerprint()`, named options forwarded to
#'   `fit_fingerprint()`, such as `max_size`, `unknown`, and `na_action`.
#'   For `fit_fingerprint()`, reserved for future options and must be empty.
#'   Options after `...` must be named.
#' @param max_size Maximum number of residues per connected substructure, an
#'   integer from 1 to 7. Includes all sizes from 1 through this maximum.
#'   Defaults to 3, including mono-, di-, and branched trisaccharides.
#' @param unknown How transform handles features absent from the fitted
#'   dictionary: `"ignore"` (default), `"warn"`, or `"error"`. Unknown features
#'   never add dimensions. Warnings are issued once per batch.
#' @param na_action `"propagate"` (default) excludes missing glycans from fitting
#'   and returns an all-`NA_integer_` vector at each missing input position.
#'   `"error"` rejects missing glycans during either fitting or transformation.
#' @param encoder An unmodified object returned by `fit_fingerprint()`.
#'
#' @details
#' Features are canonical IUPAC strings returned by [glymotif::extract_motif()].
#' Counts come directly from [glymotif::count_motifs()] with substructure
#' alignment, linkages included, strict substituent matching, and strict mode.
#' Overlapping occurrences count separately; symmetric mappings of the same
#' residue set count once. Binary values are `as.integer(count > 0)`.
#'
#' glymotif semantics apply: concrete residues match compatible generic motifs,
#' and unknown motif linkage positions act as wildcards; the reverse does not
#' hold in strict mode. Feature columns can therefore overlap semantically.
#' Motif reducing-end half-linkages constrain residue anomers, even when the
#' residue occurs internally. Linkage positions and substituents are retained.
#' Alditol structures and unresolved floating parts/substituents are rejected in
#' this version. Protein/lipid attachments, 3D geometry, and atomic stereochemistry
#' beyond the residue/linkage labels are not encoded.
#'
#' **Each top-level call fits its own batch. Independently computed batches must
#' not be compared directly, even when their vectors have equal lengths.**
#' Combine all glycans into one call, or fit once and reuse the encoder. In a
#' predictive workflow fit on training data only, then transform held-out data.
#' Unknown fragments are projected away by default; an all-zero vector can mean
#' that no known fragments occurred, not that the glycan is empty.
#'
#' The dictionary is the unique extracted motifs sorted by UTF-8 bytes, without
#' hashing or folding. Child order and input order do not affect it. Unknown
#' features mean newly extracted canonical motifs absent from this dictionary,
#' even if a broader known motif matches them. Save the encoder with `saveRDS()`
#' and restore with `readRDS()`. It stores ordinary R data and dependency versions,
#' not training graphs or external pointers. Transform requires the same
#' glyrepr, glymotif, and glyparse versions used for fitting, so changed upstream
#' normalization or matching semantics cannot silently alter saved encoders.
#'
#' An empty input returns `list()`. Fitting empty or all-missing input creates
#' an empty dictionary: each transformed vector then has length zero, including
#' missing entries, whose missingness cannot be represented at zero dimensions.
#' Keep the original `is.na(glycans)` mask when using such an encoder.
#' Large, highly branched structures can be expensive to enumerate; start with
#' the default `max_size = 3` and increase only as needed. Model embeddings are
#' not implemented in this version.
#'
#' @returns `glycan_fingerprint()` and `transform_fingerprint()` return a list
#'   of named integer vectors with the same length, order, and names as `glycans`.
#'   All vectors use `encoder$features` as their names, including zero entries.
#'   `fit_fingerprint()` returns a `glyfp_encoder` list containing `features`,
#'   `method`, `params`, `schema_version`, and `versions`.
#' @examples
#' library(glyrepr)
#' a <- as_glycan_structure(c("Gal(b1-4)GlcNAc", "GlcNAc"))
#' b <- as_glycan_structure("Gal(b1-3)GalNAc")
#'
#' # One batch provides one common feature space.
#' all_fp <- glycan_fingerprint(c(a, b), method = "substructure_binary")
#' lengths(all_fp)
#'
#' # Or retain an encoder for future batches.
#' encoder <- fit_fingerprint(a, method = "substructure_count")
#' a_fp <- transform_fingerprint(a, encoder)
#' b_fp <- transform_fingerprint(b, encoder)
#' identical(names(a_fp[[1]]), names(b_fp[[1]]))
#'
#' path <- tempfile(fileext = ".rds")
#' saveRDS(encoder, path)
#' identical(transform_fingerprint(b, readRDS(path)), b_fp)
#' unlink(path)
#' @export
glycan_fingerprint <- function(glycans, method = "substructure_binary", ...) {
  encoder <- fit_fingerprint(glycans, method, ...)
  transform_fingerprint(glycans, encoder)
}

#' @rdname glycan_fingerprint
#' @export
fit_fingerprint <- function(
  glycans,
  method = "substructure_binary",
  ...,
  max_size = 3L,
  unknown = "ignore",
  na_action = "propagate"
) {
  rlang::check_dots_empty()
  validate_options(method, max_size, unknown, na_action)
  params <- list(
    max_size = as.integer(max_size),
    unknown = unknown,
    na_action = na_action
  )
  glycans <- normalize_glycans(glycans, na_action)
  features <- extract_features(glycans, params$max_size)
  structure(
    list(
      features = features,
      method = method,
      params = params,
      schema_version = 1L,
      versions = engine_versions()
    ),
    class = "glyfp_encoder"
  )
}

#' @rdname glycan_fingerprint
#' @export
transform_fingerprint <- function(glycans, encoder) {
  validate_encoder(encoder)
  glycans <- normalize_glycans(glycans, encoder$params$na_action)
  if (!identical(encoder$versions, engine_versions())) {
    cli::cli_abort(
      "Encoder dependency versions differ from this session. Restore the recorded environment or fit a new encoder and recompute all batches.",
      class = "glyfp_error_version"
    )
  }
  if (encoder$params$unknown != "ignore") {
    unseen <- setdiff(
      extract_features(glycans, encoder$params$max_size),
      encoder$features
    )
    if (length(unseen)) {
      message <- c(
        "Found {length(unseen)} extracted motif(s) outside the encoder dictionary.",
        "i" = "No new dimensions are added; known motifs are still matched."
      )
      if (encoder$params$unknown == "error") {
        cli::cli_abort(message, class = "glyfp_error_unknown")
      }
      cli::cli_warn(message, class = "glyfp_warning_unknown")
    }
  }
  missing <- is.na(glycans)
  result <- lapply(seq_along(glycans), function(i) {
    stats::setNames(
      rep.int(if (missing[i]) NA_integer_ else 0L, length(encoder$features)),
      encoder$features
    )
  })
  valid <- which(!missing)
  if (length(valid) && length(encoder$features)) {
    # Deduplicate the batch before matching, then expand without losing names.
    unique_glycans <- unique(glycans[valid])
    motifs <- glyrepr::as_glycan_structure(encoder$features)
    counts <- glymotif::count_motifs(
      unique_glycans,
      motifs,
      alignments = "substructure",
      ignore_linkages = FALSE,
      strict_sub = TRUE,
      mode = "strict"
    )
    if (encoder$method == "substructure_binary") {
      counts <- 1L * (counts > 0L)
    }
    rows <- match(glycans[valid], unique_glycans)
    result[valid] <- lapply(rows, function(i) {
      stats::setNames(as.integer(counts[i, ]), encoder$features)
    })
  }
  names(result) <- names(glycans)
  result
}

#' @export
print.glyfp_encoder <- function(x, ...) {
  validate_encoder(x)
  cat(
    "<glyfp_encoder>\n",
    "Method: ",
    x$method,
    "\n",
    "Features: ",
    length(x$features),
    "\n",
    "Maximum residues: ",
    x$params$max_size,
    "\n",
    sep = ""
  )
  invisible(x)
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
  choice(method, c("substructure_binary", "substructure_count"), "method")
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

validate_encoder <- function(x) {
  ok <- inherits(x, "glyfp_encoder") &&
    is.list(x) &&
    identical(
      names(x),
      c("features", "method", "params", "schema_version", "versions")
    ) &&
    identical(x$schema_version, 1L) &&
    is.character(x$features) &&
    !anyNA(x$features) &&
    !anyDuplicated(x$features) &&
    all(nzchar(x$features)) &&
    identical(x$features, sort(x$features, method = "radix")) &&
    is.list(x$params) &&
    identical(names(x$params), c("max_size", "unknown", "na_action")) &&
    is.character(x$versions) &&
    identical(names(x$versions), c("glyrepr", "glymotif", "glyparse")) &&
    !anyNA(x$versions) &&
    all(nzchar(x$versions))
  if (ok) {
    ok <- tryCatch(
      {
        validate_options(
          x$method,
          x$params$max_size,
          x$params$unknown,
          x$params$na_action
        )
        TRUE
      },
      error = function(e) FALSE
    )
  }
  if (!ok) {
    cli::cli_abort(
      "{.arg encoder} must be an unmodified, supported glyfp encoder.",
      class = "glyfp_error_encoder"
    )
  }
  invisible(x)
}
