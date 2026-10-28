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
#'   `"substructure_count"` (number of occurrences), or `"sweetnet"` (optional
#'   pretrained 128-dimensional CPU embedding).
#' @param ... For `glycan_fingerprint()`, named options forwarded to
#'   `fit_fingerprint()`, such as `max_size`, `unknown`, and `na_action`.
#'   For `fit_fingerprint()`, reserved for future options and must be empty.
#'   Options after `...` must be named.
#' @param max_size Maximum number of residues per connected substructure, an
#'   integer from 1 to 7. Includes all sizes from 1 through this maximum.
#'   Defaults to 3, including mono-, di-, and branched trisaccharides. Applies
#'   only to substructure methods; supplying it for SweetNet is an error.
#' @param unknown How transform handles features absent from the fitted
#'   dictionary: `"ignore"` (default), `"warn"`, or `"error"`. Unknown features
#'   never add dimensions. Warnings are issued once per batch. SweetNet defaults
#'   to `"error"` and supports only that policy: tokens outside its fixed pretrained
#'   vocabulary are rejected, not ignored or assigned untrained representations.
#' @param na_action `"propagate"` (default) excludes missing glycans from fitting
#'   and returns an all-NA vector at each missing input position (integer for
#'   substructures, double for SweetNet).
#'   `"error"` rejects missing glycans during either fitting or transformation.
#' @param checkpoint For SweetNet, `"sweetnet-v2"` (default), `"sweetnet-v1"`
#'   (legacy), or a local path to either supported official checkpoint. Named
#'   checkpoints are downloaded from a fixed official revision on first use and
#'   verified with SHA-256 before caching or loading. Local files are identified
#'   by hash and never replaced or downloaded. Both models return 128 coordinates,
#'   but their coordinate spaces differ. This option applies only to SweetNet.
#' @param cache_dir Directory for named SweetNet checkpoints. Defaults to the
#'   per-user cache returned by `tools::R_user_dir("glyfp", "cache")`. No directory
#'   is created for a local checkpoint or an offline cache miss.
#' @param download Whether a missing named checkpoint may be downloaded (default
#'   `TRUE`). Set `FALSE` for offline use; a verified cache entry or local file is
#'   then required. Corrupt cache files are rejected, never silently overwritten.
#'   This policy and the cache directory are stored in the encoder and reused
#'   during transformation. Loading the package does not access the network.
#' @param encoder An unmodified object returned by `fit_fingerprint()`.
#'
#' @details
#' For substructure methods, features are canonical IUPAC strings returned by [glymotif::extract_motif()].
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
#' **Each top-level substructure call fits its own batch. Independently computed batches must
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
#' An empty input returns `list()`. For substructure methods, fitting empty or all-missing input creates
#' an empty dictionary: each transformed vector then has length zero, including
#' missing entries, whose missingness cannot be represented at zero dimensions.
#' Keep the original `is.na(glycans)` mask when using such an encoder.
#' Large, highly branched structures can be expensive to enumerate; start with
#' the default `max_size = 3` and increase only as needed.
#'
#' SweetNet supports the legacy and v2 pretrained species models from glycowork.
#' The default is v2, matching the current upstream pretrained model. Fitting binds
#' a fixed pretrained coordinate space and validates inputs, without training.
#' It always has 128 features, including when fitted to empty/all-missing input.
#' Independent SweetNet fits are comparable only when their recorded model and
#' preprocessing metadata match. Reusing one encoder is recommended.
#' Unknown reducing-end anomers at the residue's inferred normal carbon position
#' are omitted from the model input; this is glyrepr's default annotation for an
#' unspecified reducing end. Explicit alpha/beta half-linkages, other reducing-end
#' positions, and ambiguous internal linkages are rejected. The directed
#' residue/linkage graph is built directly from the glyrepr structure.
#' No unknown model token is accepted. See `vignette("glyfp", package = "glyfp")`.
#'
#' SweetNet runs entirely in R using the optional torch, digest, and zip packages.
#' Its encoder records model version, official source revision, weight and
#' vocabulary SHA-256, preprocessing/readout, R implementation hashes and the R
#' torch version. Every transform validates these records. Independently fitted
#' v1 and v2 encoders must not be mixed, even though both have 128 features with
#' the same coordinate names. Retain the encoder alongside exported vectors.
#' Named checkpoints can be reacquired if the cache is removed and downloading
#' is allowed; local checkpoints must remain at their recorded paths. No model
#' pointers or weights are serialized. The legacy archive stores CUDA tags;
#' glyfp relocates these tags to CPU in a temporary copy before loading. Tensor
#' storage bytes and the original file are unchanged. No Python is used.
#' The optional torch runtime must be installed explicitly by the user; it is
#' never automatically installed by these functions.
#' CPU evaluation uses float32 tensors, one graph at a time, frozen normalization,
#' no dropout and no gradients, without changing the caller's thread settings.
#' Cross-platform/runtime comparisons should allow floating-point tolerance
#' (tested at a maximum absolute coordinate error of 1e-5).
#'
#' @returns `glycan_fingerprint()` and `transform_fingerprint()` return a list
#'   of named integer vectors (substructures) or double vectors (SweetNet), with
#'   the same length, order, and names as `glycans`. Missing SweetNet inputs yield
#'   128 `NA_real_` values.
#'   All vectors use `encoder$features` as their names, including zero entries.
#'   `fit_fingerprint()` returns a `glyfp_encoder` list containing `features`,
#'   `method`, `params`, `schema_version`, and `versions`. For SweetNet, `features`
#'   are `sweetnet_001` through `sweetnet_128`; `params$model` fixes model metadata.
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
  na_action = "propagate",
  checkpoint = "sweetnet-v2",
  cache_dir = tools::R_user_dir("glyfp", "cache"),
  download = TRUE
) {
  rlang::check_dots_empty()
  if (identical(method, "sweetnet") && missing(unknown)) {
    unknown <- "error"
  }
  if (identical(method, "sweetnet") && !missing(max_size)) {
    cli::cli_abort(
      "{.arg max_size} applies only to substructure methods.",
      class = "glyfp_error_argument"
    )
  }
  if (
    !identical(method, "sweetnet") &&
      (!missing(checkpoint) || !missing(cache_dir) || !missing(download))
  ) {
    cli::cli_abort(
      "{.arg checkpoint}, {.arg cache_dir}, and {.arg download} apply only to SweetNet.",
      class = "glyfp_error_argument"
    )
  }
  validate_options(method, max_size, unknown, na_action)
  params <- list(
    max_size = as.integer(max_size),
    unknown = unknown,
    na_action = na_action
  )
  glycans <- normalize_glycans(glycans, na_action)
  if (identical(method, "sweetnet")) {
    return(fit_sweetnet(
      glycans,
      checkpoint,
      unknown,
      na_action,
      cache_dir,
      download
    ))
  }
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
  if (identical(encoder$method, "sweetnet")) {
    return(transform_sweetnet(glycans, encoder))
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
    sep = ""
  )
  if (identical(x$method, "sweetnet")) {
    cat("Model: ", x$params$model$model, " (CPU, 128 dimensions)\n", sep = "")
  } else {
    cat("Maximum residues: ", x$params$max_size, "\n", sep = "")
  }
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
  choice(
    method,
    c("substructure_binary", "substructure_count", "sweetnet"),
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
    (if (identical(x$method, "sweetnet")) {
      valid_sweetnet_params(x$params) &&
        identical(x$features, sweetnet_features())
    } else {
      identical(names(x$params), c("max_size", "unknown", "na_action"))
    }) &&
    is.character(x$versions) &&
    identical(names(x$versions), c("glyrepr", "glymotif", "glyparse")) &&
    !anyNA(x$versions) &&
    all(nzchar(x$versions))
  if (ok) {
    ok <- tryCatch(
      {
        validate_options(
          x$method,
          if (identical(x$method, "sweetnet")) 3L else x$params$max_size,
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
