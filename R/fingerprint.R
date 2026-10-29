validate_choice <- function(x, choices, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !x %in% choices) {
    cli::cli_abort(
      "{.arg {arg}} must be one of {.val {choices}}.",
      class = "glyfp_error_argument"
    )
  }
}

complete_glycans <- function(glycans, na_action) {
  validate_choice(na_action, c("propagate", "error"), "na_action")
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
  glycans
}

embedding_fingerprints <- function(glycans, run, features, model) {
  glycans <- embedding_projection(glycans)
  valid <- which(!is.na(glycans))
  unique_glycans <- unique(glycans[valid])
  values <- run(unique_glycans)
  if (
    length(values) != length(unique_glycans) ||
      any(vapply(
        values,
        \(x) length(x) != length(features) || any(!is.finite(x)),
        TRUE
      ))
  ) {
    cli::cli_abort(
      "{model} returned invalid embeddings.",
      class = "glyfp_error_backend"
    )
  }
  result <- lapply(seq_along(glycans), \(i) {
    stats::setNames(rep(NA_real_, length(features)), features)
  })
  result[valid] <- lapply(match(glycans[valid], unique_glycans), \(i) {
    stats::setNames(values[[i]], features)
  })
  names(result) <- names(glycans)
  result
}
