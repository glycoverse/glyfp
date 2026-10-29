# Mirror glymotif's public input contract without calling its private helpers.
normalize_glycans <- function(glycans, na_action) {
  original_names <- names(glycans)
  if (!inherits(glycans, "glyrepr_structure")) {
    if (!is.character(glycans) || !is.null(dim(glycans))) {
      cli::cli_abort(
        "{.arg glycans} must be a glyrepr structure vector or glycan structure character vector.",
        class = "glyfp_error_input"
      )
    }
    glycans <- tryCatch(
      glyparse::auto_parse(glycans, on_failure = "error"),
      error = function(e) {
        cli::cli_abort(
          "Some {.arg glycans} could not be parsed as valid glycan structures.",
          parent = e,
          class = "glyfp_error_parse"
        )
      }
    )
  }
  if (na_action == "error" && anyNA(glycans)) {
    cli::cli_abort(
      "Missing glycans at position(s): {paste(which(is.na(glycans)), collapse = ', ')}.",
      class = "glyfp_error_missing"
    )
  }
  # Validate through glyrepr's public constructor, including corrupted storage.
  graphs <- tryCatch(
    as.list(unique(glycans[!is.na(glycans)])),
    error = function(e) {
      cli::cli_abort(
        "Invalid glyrepr structure storage.",
        parent = e,
        class = "glyfp_error_input"
      )
    }
  )
  for (g in graphs) {
    if (is.null(g)) {
      cli::cli_abort(
        "Missing graph for a non-missing glycan.",
        class = "glyfp_error_input"
      )
    }
    tryCatch(glyrepr::glycan_structure(g), error = function(e) {
      cli::cli_abort(
        "Invalid glycan structure.",
        parent = e,
        class = "glyfp_error_input"
      )
    })
    if (
      length(igraph::graph_attr(g, "floating_parts")) ||
        length(igraph::graph_attr(g, "floating_substituents")) ||
        isTRUE(igraph::graph_attr(g, "alditol"))
    ) {
      cli::cli_abort(
        "Unresolved floating parts/substituents and alditols are not supported in this version.",
        class = "glyfp_error_unsupported"
      )
    }
  }
  names(glycans) <- original_names
  glycans
}

extract_features <- function(glycans, max_size) {
  glycans <- unique(glycans[!is.na(glycans)])
  if (!length(glycans)) {
    return(character())
  }
  motifs <- glymotif::extract_motif(glycans, max_size = max_size)
  sort(unique(enc2utf8(as.character(motifs))), method = "radix")
}

.motif_dictionary <- new.env(parent = emptyenv())

motif_features <- function(max_size) {
  if (is.null(.motif_dictionary$data)) {
    .motif_dictionary$data <- utils::read.delim(
      system.file(
        "extdata",
        "motif-dictionary.tsv",
        package = "glyfp",
        mustWork = TRUE
      ),
      quote = "",
      comment.char = "",
      colClasses = c("integer", "character")
    )
  }
  dictionary <- .motif_dictionary$data
  dictionary$motif[dictionary$size <= max_size]
}

motif_fingerprints <- function(glycans, method, max_size, unknown) {
  features <- motif_features(max_size)
  valid <- which(!is.na(glycans))
  unique_glycans <- unique(glycans[valid])
  extracted <- extract_features(unique_glycans, max_size)
  unseen <- setdiff(extracted, features)
  if (length(unseen) && unknown != "ignore") {
    message <- c(
      "Found {length(unseen)} motif(s) outside the built-in glydb dictionary.",
      "i" = "No dimensions are added; known motifs are still counted."
    )
    if (unknown == "error") {
      cli::cli_abort(message, class = "glyfp_error_unknown")
    }
    cli::cli_warn(message, class = "glyfp_warning_unknown")
  }
  # With concrete, fully specified motifs, any matching connected subgraph is
  # present in this enumeration. Match only those candidates, then fill the
  # fixed dictionary with zeros. This preserves overlapping glymotif counts.
  selected <- which(features %in% extracted)
  counts <- matrix(0L, nrow = length(unique_glycans), ncol = length(features))
  if (length(selected) && length(unique_glycans)) {
    counts[, selected] <- glymotif::count_motifs(
      unique_glycans,
      glyrepr::as_glycan_structure(features[selected]),
      alignments = "substructure",
      ignore_linkages = FALSE,
      strict_sub = TRUE,
      mode = "strict"
    )
  }
  if (method == "substructure_binary") {
    counts <- 1L * (counts > 0L)
  }
  result <- lapply(seq_along(glycans), function(i) {
    stats::setNames(rep(NA_integer_, length(features)), features)
  })
  result[valid] <- lapply(match(glycans[valid], unique_glycans), function(i) {
    stats::setNames(as.integer(counts[i, ]), features)
  })
  names(result) <- names(glycans)
  result
}
