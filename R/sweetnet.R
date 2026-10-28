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
    # Build the directed bipartite graph directly from glyrepr, preserving
    # residue/substituent labels and each parent -> linkage -> child edge.
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

sweetnet_path <- function(value, arg) {
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
  if (!file.exists(value) || dir.exists(value)) {
    cli::cli_abort(
      "{.arg {arg}} does not identify an existing local file.",
      class = "glyfp_error_backend"
    )
  }
  normalizePath(value, mustWork = TRUE)
}

sweetnet_require_namespace <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

sweetnet_restore_rng <- function() {
  seeded <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (seeded) {
    seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  }
  function() {
    if (seeded) {
      assign(".Random.seed", seed, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }
}

sweetnet_dependencies <- function() {
  # torch's namespace startup can auto-install its runtime and seed its generator
  # from R's RNG. Loading this optional backend must do neither on the caller's behalf.
  install <- Sys.getenv("TORCH_INSTALL", unset = NA_character_)
  restore_rng <- sweetnet_restore_rng()
  on.exit(restore_rng(), add = TRUE)
  on.exit(
    {
      if (is.na(install)) {
        Sys.unsetenv("TORCH_INSTALL")
      } else {
        Sys.setenv(TORCH_INSTALL = install)
      }
    },
    add = TRUE
  )
  Sys.setenv(TORCH_INSTALL = "0")
  for (pkg in c("torch", "digest", "zip")) {
    if (!sweetnet_require_namespace(pkg)) {
      cli::cli_abort(
        "Install {.pkg {pkg}} to use SweetNet.",
        class = "glyfp_error_backend"
      )
    }
  }
  if (!torch::torch_is_installed()) {
    cli::cli_abort(
      "Install the R torch runtime with {.code torch::install_torch()} to use SweetNet.",
      class = "glyfp_error_backend"
    )
  }
}

sweetnet_sha256 <- function(path) digest::digest(file = path, algo = "sha256")

sweetnet_vocabulary <- function() {
  path <- system.file(
    "extdata",
    "sweetnet-vocabulary.tsv",
    package = "glyfp",
    mustWork = TRUE
  )
  if (
    !identical(
      sweetnet_sha256(path),
      "9edf1dbdb4df99ae6444108ba09e72f786c175c0b051a32fd82af3e46edc69fc"
    )
  ) {
    cli::cli_abort(
      "SweetNet vocabulary resource has changed.",
      class = "glyfp_error_backend"
    )
  }
  x <- utils::read.delim(
    path,
    header = FALSE,
    quote = "",
    comment.char = "",
    colClasses = c("character", "integer")
  )
  stats::setNames(x[[2L]], x[[1L]])
}

# This is a relocation of one hash-verified artifact, not a general pickle reader.
# Torch for R currently has no map_location option on load_state_dict(). Only
# BINUNICODE CUDA storage tags change; tensor storage bytes do not.
sweetnet_load_weights <- function(checkpoint, id) {
  spec <- sweetnet_checkpoints()[[id]]
  if (
    !identical(
      sweetnet_sha256(checkpoint),
      spec$sha256
    )
  ) {
    cli::cli_abort(
      "Checkpoint SHA-256 does not match the supported official SweetNet artifact.",
      class = "glyfp_error_backend"
    )
  }
  directory <- tempfile("glyfp-sweetnet-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  utils::unzip(checkpoint, exdir = directory)
  records <- utils::unzip(checkpoint, list = TRUE)$Name
  record <- records[grepl("^[^/]+/data[.]pkl$", records)]
  if (length(record) != 1L) {
    cli::cli_abort(
      "Unexpected checkpoint archive metadata.",
      class = "glyfp_error_backend"
    )
  }
  archive <- dirname(record)
  pickle <- file.path(directory, record)
  bytes <- readBin(pickle, "raw", n = file.info(pickle)$size)
  needle <- c(as.raw(c(88, 6, 0, 0, 0)), charToRaw("cuda:0"))
  positions <- which(vapply(
    seq_len(length(bytes) - length(needle) + 1L),
    function(i) {
      identical(bytes[i + seq_along(needle) - 1L], needle)
    },
    TRUE
  ))
  if (!length(positions)) {
    return(torch::load_state_dict(checkpoint))
  }
  for (position in rev(positions)) {
    bytes <- c(
      bytes[seq_len(position - 1L)],
      as.raw(c(88, 3, 0, 0, 0)),
      charToRaw("cpu"),
      bytes[seq.int(position + length(needle), length(bytes))]
    )
  }
  writeBin(bytes, pickle)
  relocated <- file.path(directory, "cpu.pt")
  zip::zipr(relocated, archive, root = directory)
  torch::load_state_dict(relocated)
}

sweetnet_metadata <- function(id) {
  spec <- sweetnet_checkpoints()[[id]]
  sources <- list(
    checkpoints = sweetnet_checkpoints,
    dependencies = sweetnet_dependencies,
    rng = sweetnet_restore_rng,
    preprocessing = sweetnet_inputs,
    vocabulary = sweetnet_vocabulary,
    weights = sweetnet_load_weights,
    inference = sweetnet_infer,
    execution = sweetnet_run
  )
  list(
    model = id,
    source_revision = spec$revision,
    checkpoint_sha256 = spec$sha256,
    vocabulary_sha256 = "838281860e3bbddc0de91ac635bf1ba746d79bf62891e5b86ca08eff537a974e",
    preprocessing = "glyfp-sweetnet-r-v2:unknown-inferred-root-only;exact-internal-links;directed-residue-linkage-graph",
    readout = "bn2(lin2):before-activation-dropout-classifier",
    dimensions = 128L,
    device = "cpu",
    batch_size = 1L,
    versions = list(torch = as.character(utils::packageVersion("torch"))),
    code_sha256 = lapply(sources, function(f) {
      digest::digest(
        paste(deparse(body(f)), collapse = "\n"),
        algo = "sha256",
        serialize = FALSE
      )
    })
  )
}

sweetnet_infer <- function(input, weights, vocabulary, id) {
  tokens <- unname(vocabulary[as.character(input$labels)])
  if (anyNA(tokens)) {
    unknown <- unique(as.character(input$labels)[is.na(tokens)])
    cli::cli_abort(
      "Unsupported pretrained vocabulary token(s): {.val {unknown}}.",
      class = "glyfp_error_unknown"
    )
  }
  torch::with_no_grad({
    # R torch tensor indices are one-based, while the original vocabulary is zero-based.
    indices <- torch::torch_tensor(
      tokens + 1L,
      dtype = torch::torch_long(),
      device = "cpu"
    )
    x <- torch::nnf_embedding(indices, weights[["item_embedding.weight"]])
    for (layer in seq_len(3L)) {
      aggregated <- torch::torch_zeros_like(x)
      if (length(input$edges)) {
        edges <- do.call(rbind, input$edges) + 1L
        source <- torch::torch_tensor(
          edges[, 1L],
          dtype = torch::torch_long(),
          device = "cpu"
        )
        target <- torch::torch_tensor(
          edges[, 2L],
          dtype = torch::torch_long(),
          device = "cpu"
        )
        aggregated$index_add_(1L, target, x$index_select(1L, source))
      }
      prefix <- paste0("conv", layer)
      previous <- x
      x <- torch::nnf_linear(
        aggregated,
        weights[[paste0(prefix, ".lin_rel.weight")]],
        weights[[paste0(prefix, ".lin_rel.bias")]]
      ) +
        torch::nnf_linear(x, weights[[paste0(prefix, ".lin_root.weight")]])
      if (id == "sweetnet-v2") {
        norm <- paste0("norms.", layer - 1L)
        x <- torch::nnf_batch_norm(
          x,
          weights[[paste0(norm, ".running_mean")]],
          weights[[paste0(norm, ".running_var")]],
          weights[[paste0(norm, ".weight")]],
          weights[[paste0(norm, ".bias")]],
          training = FALSE,
          eps = 1e-5
        )
      }
      x <- torch::nnf_leaky_relu(x, negative_slope = 0.01)
      if (id == "sweetnet-v2") x <- previous + x
    }
    x <- if (id == "sweetnet-v2") {
      torch::torch_cat(
        list(
          x$mean(dim = 1L, keepdim = TRUE),
          x$amax(dim = 1L, keepdim = TRUE)
        ),
        dim = 2L
      )
    } else {
      x$mean(dim = 1L, keepdim = TRUE)
    }
    for (layer in seq_len(2L)) {
      x <- torch::nnf_linear(
        x,
        weights[[paste0("lin", layer, ".weight")]],
        weights[[paste0("lin", layer, ".bias")]]
      )
      prefix <- paste0("bn", layer)
      x <- torch::nnf_batch_norm(
        x,
        weights[[paste0(prefix, ".running_mean")]],
        weights[[paste0(prefix, ".running_var")]],
        weights[[paste0(prefix, ".weight")]],
        weights[[paste0(prefix, ".bias")]],
        training = FALSE,
        eps = 1e-5
      )
      if (layer == 1L) x <- torch::nnf_leaky_relu(x, negative_slope = 0.01)
    }
    as.double(torch::as_array(x))
  })
}

sweetnet_run <- function(
  inputs,
  checkpoint,
  metadata = NULL,
  fit = FALSE,
  cache_dir = tools::R_user_dir("glyfp", "cache"),
  download = TRUE
) {
  restore_rng <- sweetnet_restore_rng()
  on.exit(restore_rng(), add = TRUE)
  sweetnet_dependencies()
  if (
    !is.null(metadata) &&
      !identical(metadata, sweetnet_metadata(metadata$model))
  ) {
    cli::cli_abort(
      "Encoder model, preprocessing, code, or dependency metadata differs; restore the fitted environment.",
      class = "glyfp_error_backend"
    )
  }
  vocabulary <- sweetnet_vocabulary()
  for (input in inputs) {
    unknown <- setdiff(as.character(input$labels), names(vocabulary))
    if (length(unknown)) {
      cli::cli_abort(
        "Unsupported pretrained vocabulary token(s): {.val {unknown}}.",
        class = "glyfp_error_unknown"
      )
    }
  }
  resolved <- sweetnet_resolve(checkpoint, cache_dir, download)
  current <- sweetnet_metadata(resolved$id)
  if (!is.null(metadata) && !identical(metadata, current)) {
    cli::cli_abort(
      "Encoder model, preprocessing, code, or dependency metadata differs; restore the fitted environment.",
      class = "glyfp_error_backend"
    )
  }
  weights <- sweetnet_load_weights(resolved$path, resolved$id)
  embeddings <- if (fit) {
    list()
  } else {
    lapply(
      inputs,
      sweetnet_infer,
      weights = weights,
      vocabulary = vocabulary,
      id = resolved$id
    )
  }
  list(
    metadata = current,
    embeddings = embeddings,
    checkpoint = resolved$source
  )
}

fit_sweetnet <- function(
  glycans,
  checkpoint,
  unknown,
  na_action,
  cache_dir,
  download
) {
  if (!identical(unknown, "error")) {
    cli::cli_abort(
      "SweetNet supports only {.code unknown = 'error'}; untrained/unknown tokens cannot be ignored.",
      class = "glyfp_error_argument"
    )
  }
  if (
    !is.character(checkpoint) ||
      length(checkpoint) != 1L ||
      is.na(checkpoint) ||
      !nzchar(checkpoint)
  ) {
    cli::cli_abort(
      "{.arg checkpoint} must be a registered name or local path.",
      class = "glyfp_error_argument"
    )
  }
  if (!is.logical(download) || length(download) != 1L || is.na(download)) {
    cli::cli_abort(
      "{.arg download} must be TRUE or FALSE.",
      class = "glyfp_error_argument"
    )
  }
  cache_dir <- sweetnet_cache_dir(cache_dir)
  inputs <- sweetnet_inputs(unique(glycans[!is.na(glycans)]))
  response <- sweetnet_run(
    inputs,
    checkpoint,
    fit = TRUE,
    cache_dir = cache_dir,
    download = download
  )
  structure(
    list(
      features = sweetnet_features(),
      method = "sweetnet",
      params = list(
        unknown = unknown,
        na_action = na_action,
        checkpoint = response$checkpoint,
        cache_dir = sweetnet_cache_dir(cache_dir),
        download = download,
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
    encoder$params$checkpoint,
    metadata = encoder$params$model,
    cache_dir = encoder$params$cache_dir,
    download = encoder$params$download
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
      c("unknown", "na_action", "checkpoint", "cache_dir", "download", "model")
    ) &&
    identical(x$unknown, "error") &&
    is.character(x$na_action) &&
    length(x$na_action) == 1L &&
    !is.na(x$na_action) &&
    x$na_action %in% c("propagate", "error") &&
    is.logical(x$download) &&
    length(x$download) == 1L &&
    !is.na(x$download) &&
    all(vapply(
      x[c("checkpoint", "cache_dir")],
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
        "source_revision",
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
    is.character(x$model$model) &&
    length(x$model$model) == 1L &&
    !is.na(x$model$model) &&
    x$model$model %in% names(sweetnet_checkpoints()) &&
    (!x$checkpoint %in% names(sweetnet_checkpoints()) ||
      identical(x$checkpoint, x$model$model)) &&
    identical(x$model$dimensions, 128L) &&
    is.list(x$model$versions) &&
    is.list(x$model$code_sha256)
}
