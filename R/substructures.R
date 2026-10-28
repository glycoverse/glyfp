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

engine_versions <- function() {
  vapply(
    c("glyrepr", "glymotif", "glyparse"),
    function(p) as.character(utils::packageVersion(p)),
    ""
  )
}
