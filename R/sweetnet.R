sweetnet_features <- function() sprintf("sweetnet_%03d", seq_len(128L))

sweetnet_inputs <- function(glycans) {
  unname(lapply(as.list(glycans), function(g) {
    root <- which(igraph::degree(g, mode = "in") == 0L)
    mono <- igraph::vertex_attr(g, "mono")
    anomer <- igraph::graph_attr(g, "anomer")
    inferred <- paste0("?", glyrepr::infer_anomer_pos(mono[root]))
    if (length(root) != 1L || !identical(anomer, inferred)) {
      cli::cli_abort(
        "SweetNet requires an unspecified reducing-end anomer with its normal inferred carbon position; explicit alpha/beta or other half-linkages cannot be discarded.",
        class = "glyfp_error_unsupported"
      )
    }
    linkage <- igraph::edge_attr(g, "linkage")
    if (any(!grepl("^[ab][0-9]+-[0-9]+$", linkage))) {
      cli::cli_abort(
        "SweetNet currently requires fully specified internal linkages.",
        class = "glyfp_error_unsupported"
      )
    }
    iupac <- as.character(glyrepr::glycan_structure(g))
    suffix <- paste0("(", anomer, "-")
    if (!endsWith(iupac, suffix)) {
      cli::cli_abort(
        "Unexpected reducing-end serialization.",
        class = "glyfp_error_input"
      )
    }
    iupac <- substr(iupac, 1L, nchar(iupac) - nchar(suffix))
    substituent <- igraph::vertex_attr(g, "sub")
    substituent[is.na(substituent)] <- ""
    labels <- c(paste0(mono, gsub(",", "", substituent, fixed = TRUE)), linkage)
    # Compare the model parser's labeled, directed residue/linkage graph
    # against the original glyrepr graph before computing any representation.
    endpoints <- igraph::as_edgelist(g, names = FALSE)
    n <- length(mono)
    edges <- unlist(
      lapply(seq_along(linkage), function(i) {
        list(
          c(endpoints[i, 1L] - 1L, n + i - 1L),
          c(n + i - 1L, endpoints[i, 2L] - 1L)
        )
      }),
      recursive = FALSE
    )
    if (is.null(edges)) {
      edges <- list()
    }
    list(iupac = iupac, labels = I(labels), edges = edges)
  }))
}

sweetnet_path <- function(value, arg, executable = FALSE) {
  if (
    !is.character(value) ||
      length(value) != 1L ||
      is.na(value) ||
      !nzchar(value)
  ) {
    cli::cli_abort(
      "{.arg {arg}} must name an existing local file.",
      class = "glyfp_error_argument"
    )
  }
  if (executable && !file.exists(value)) {
    value <- unname(Sys.which(value))
  }
  if (!nzchar(value) || !file.exists(value) || dir.exists(value)) {
    cli::cli_abort(
      "{.arg {arg}} does not identify an existing local file.",
      class = "glyfp_error_backend"
    )
  }
  normalizePath(value, mustWork = TRUE)
}

sweetnet_run <- function(
  inputs,
  python,
  checkpoint,
  metadata = NULL,
  fit = FALSE
) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    cli::cli_abort(
      "Install {.pkg jsonlite} to use SweetNet.",
      class = "glyfp_error_backend"
    )
  }
  bridge <- system.file(
    "python",
    "sweetnet.py",
    package = "glyfp",
    mustWork = TRUE
  )
  paths <- vapply(seq_len(4L), function(i) tempfile(), "")
  on.exit(unlink(paths), add = TRUE)
  jsonlite::write_json(
    list(
      glycans = inputs,
      checkpoint = checkpoint,
      expected_metadata = metadata,
      fit = fit,
      preprocessing_source = paste(
        deparse(body(sweetnet_inputs)),
        collapse = "\n"
      )
    ),
    paths[1L],
    auto_unbox = TRUE,
    null = "null",
    digits = NA
  )
  status <- tryCatch(
    suppressWarnings(system2(
      python,
      c(shQuote(bridge), shQuote(paths[1L]), shQuote(paths[2L])),
      stdout = paths[3L],
      stderr = paths[4L]
    )),
    error = function(e) {
      cli::cli_abort(
        "Could not start SweetNet Python.",
        parent = e,
        class = "glyfp_error_backend"
      )
    }
  )
  response <- if (file.exists(paths[2L])) {
    tryCatch(
      jsonlite::read_json(paths[2L], simplifyVector = FALSE),
      error = function(e) NULL
    )
  } else {
    NULL
  }
  if (!identical(status, 0L) || is.null(response) || !isTRUE(response$ok)) {
    message <- if (
      is.character(response$error) && length(response$error) == 1L
    ) {
      response$error
    } else {
      paste(
        utils::tail(readLines(paths[4L], warn = FALSE), 8L),
        collapse = "\n"
      )
    }
    if (!nzchar(message)) {
      message <- "The Python process did not return a valid result."
    }
    cli::cli_abort(
      "SweetNet backend: {message}",
      class = if (identical(response$error_type, "unknown")) {
        "glyfp_error_unknown"
      } else {
        "glyfp_error_backend"
      }
    )
  }
  response
}

fit_sweetnet <- function(glycans, python, checkpoint, unknown, na_action) {
  if (!identical(unknown, "error")) {
    cli::cli_abort(
      "SweetNet supports only {.code unknown = 'error'}; untrained/unknown tokens cannot be ignored.",
      class = "glyfp_error_argument"
    )
  }
  python <- sweetnet_path(python, "python", executable = TRUE)
  checkpoint <- sweetnet_path(checkpoint, "checkpoint")
  inputs <- sweetnet_inputs(unique(glycans[!is.na(glycans)]))
  response <- sweetnet_run(inputs, python, checkpoint, fit = TRUE)
  structure(
    list(
      features = sweetnet_features(),
      method = "sweetnet",
      params = list(
        unknown = unknown,
        na_action = na_action,
        python = python,
        checkpoint = checkpoint,
        model = response$metadata
      ),
      schema_version = 1L,
      versions = engine_versions()
    ),
    class = "glyfp_encoder"
  )
}

transform_sweetnet <- function(glycans, encoder) {
  missing <- is.na(glycans)
  valid <- which(!missing)
  result <- lapply(seq_along(glycans), function(i) {
    stats::setNames(rep(NA_real_, 128L), encoder$features)
  })
  unique_glycans <- unique(glycans[valid])
  inputs <- sweetnet_inputs(unique_glycans)
  # Even empty batches validate the fixed environment and artifact.
  response <- sweetnet_run(
    inputs,
    encoder$params$python,
    encoder$params$checkpoint,
    metadata = encoder$params$model
  )
  values <- lapply(response$embeddings, function(row) as.double(unlist(row)))
  if (
    length(values) != length(unique_glycans) ||
      any(vapply(
        values,
        function(row) length(row) != 128L || any(!is.finite(row)),
        TRUE
      ))
  ) {
    cli::cli_abort(
      "SweetNet returned invalid embedding dimensions or values.",
      class = "glyfp_error_backend"
    )
  }
  result[valid] <- lapply(match(glycans[valid], unique_glycans), function(i) {
    stats::setNames(values[[i]], encoder$features)
  })
  names(result) <- names(glycans)
  result
}

valid_sweetnet_params <- function(x) {
  is.list(x) &&
    identical(
      names(x),
      c("unknown", "na_action", "python", "checkpoint", "model")
    ) &&
    identical(x$unknown, "error") &&
    is.character(x$na_action) &&
    length(x$na_action) == 1L &&
    !is.na(x$na_action) &&
    x$na_action %in% c("propagate", "error") &&
    all(vapply(
      x[c("python", "checkpoint")],
      function(p) {
        is.character(p) && length(p) == 1L && !is.na(p) && nzchar(p)
      },
      TRUE
    )) &&
    is.list(x$model) &&
    identical(
      names(x$model),
      c(
        "model",
        "checkpoint_sha256",
        "vocabulary_sha256",
        "preprocessing",
        "readout",
        "dimensions",
        "device",
        "batch_size",
        "versions",
        "code_sha256"
      )
    ) &&
    identical(x$model$model, "sweetnet-glycowork-1.9.0-species") &&
    identical(x$model$dimensions, 128L) &&
    is.list(x$model$versions) &&
    is.list(x$model$code_sha256)
}
