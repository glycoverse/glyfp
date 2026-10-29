# Official GlycanGT large / node-and-edge masking 35% checkpoint. Local only.
glycangt_spec <- function() {
  list(
    model = "glycangt-large-mask35",
    source_revision = "96611518c971deb89215ca163deaf9de3a59fa32",
    checkpoint_revision = "7c9f4a19b366ca48ed03995a0ad47fa1e745f24c",
    sha256 = "d490d18795ee6d4746b9b4abd0a5c733d8f7e4261998d050c1b741b65e530a54",
    vocabulary_sha256 = c(
      monomer = "3e6053ef469d957f07e3f39374357fb6f328827d8683908c8d02f54ab561a644",
      linkage = "a8e103bd141fd71825ba58749b38fef1d263fe68d1f01a6bea35665dd969ab1f"
    ),
    config = list(
      dimensions = 768L,
      heads = 32L,
      layers = 12L,
      orf_dimensions = 64L
    )
  )
}

glycangt_features <- function() sprintf("glycangt_%03d", seq_len(768L))

glycangt_checkpoint <- function(checkpoint) {
  hint <- "Download the official GlycanGT large mask-35 checkpoint yourself from https://huggingface.co/Akikitani295/GlycanGT and pass its local path as `checkpoint`. glyfp never downloads GlycanGT weights."
  if (
    !is.character(checkpoint) ||
      length(checkpoint) != 1L ||
      is.na(checkpoint) ||
      !nzchar(checkpoint) ||
      grepl("^[A-Za-z][A-Za-z0-9+.-]*://", checkpoint)
  ) {
    cli::cli_abort(
      c("GlycanGT requires a local {.arg checkpoint} file.", "i" = hint),
      class = "glyfp_error_argument"
    )
  }
  if (!file.exists(checkpoint) || dir.exists(checkpoint)) {
    cli::cli_abort(
      c(
        "GlycanGT checkpoint file does not exist: {.path {checkpoint}}.",
        "i" = hint
      ),
      class = "glyfp_error_backend"
    )
  }
  normalizePath(checkpoint, mustWork = TRUE)
}

glycangt_dependencies <- function() {
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
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cli::cli_abort(
        "Install {.pkg {pkg}} to use GlycanGT.",
        class = "glyfp_error_backend"
      )
    }
  }
  if (!torch::torch_is_installed()) {
    cli::cli_abort(
      "Install the R torch runtime with {.code torch::install_torch()} to use GlycanGT.",
      class = "glyfp_error_backend"
    )
  }
}

glycangt_vocabulary <- function() {
  hashes <- glycangt_spec()$vocabulary_sha256
  lapply(stats::setNames(names(hashes), names(hashes)), function(name) {
    path <- system.file(
      "extdata",
      paste0("glycangt-", name, ".tsv"),
      package = "glyfp",
      mustWork = TRUE
    )
    if (!identical(sweetnet_sha256(path), unname(hashes[name]))) {
      cli::cli_abort(
        "GlycanGT vocabulary resource has changed.",
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
  })
}

# Cache only tensors in memory. Every use still hashes the local original file.
# Neither checkpoints nor credentials are downloaded or persisted by this backend.
.glycangt_state <- new.env(parent = emptyenv())

glycangt_load_weights <- function(checkpoint) {
  if (!identical(sweetnet_sha256(checkpoint), glycangt_spec()$sha256)) {
    cli::cli_abort(
      "Checkpoint SHA-256 does not match the supported official GlycanGT large mask-35 artifact.",
      class = "glyfp_error_backend"
    )
  }
  if (!is.null(.glycangt_state$weights)) {
    return(.glycangt_state$weights)
  }
  directory <- tempfile("glyfp-glycangt-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  records <- utils::unzip(checkpoint, list = TRUE)$Name
  record <- records[grepl("^[^/]+/data[.]pkl$", records)]
  if (length(record) != 1L) {
    cli::cli_abort(
      "Unexpected GlycanGT archive metadata.",
      class = "glyfp_error_backend"
    )
  }
  utils::unzip(checkpoint, exdir = directory)
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
    cli::cli_abort(
      "Unexpected GlycanGT device metadata.",
      class = "glyfp_error_backend"
    )
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
  zip::zipr(
    relocated,
    dirname(record),
    root = directory,
    compression_level = 0L
  )
  weights <- torch::load_state_dict(relocated)
  glycangt_validate_weights(weights, glycangt_spec()$config)
  .glycangt_state$weights <- weights
  weights
}

glycangt_metadata <- function() {
  spec <- glycangt_spec()
  sources <- list(
    spec = glycangt_spec,
    dependencies = glycangt_dependencies,
    rng = sweetnet_restore_rng,
    preprocessing = glycangt_inputs,
    vocabulary = glycangt_vocabulary,
    checkpoint = glycangt_checkpoint,
    weights = glycangt_load_weights,
    validation = glycangt_validate_weights,
    policy = glycangt_policy,
    policy_validation = glycangt_validate_policy,
    graph_seed = glycangt_graph_seed,
    identifiers = glycangt_orf,
    inference = glycangt_infer,
    execution = glycangt_run
  )
  list(
    model = spec$model,
    source_revision = spec$source_revision,
    checkpoint_revision = spec$checkpoint_revision,
    checkpoint_sha256 = spec$sha256,
    vocabulary_sha256 = spec$vocabulary_sha256,
    config = spec$config,
    preprocessing = "glyrepr-canonical:concrete-linkages:unspecified-root:reverse-triples",
    readout = "final-graph-token",
    policy = glycangt_policy(),
    device = "cpu",
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

glycangt_backend <- function(
  glycans,
  checkpoint,
  metadata = NULL,
  fit = FALSE
) {
  # Local path validation precedes any optional dependency startup.
  checkpoint <- glycangt_checkpoint(checkpoint)
  restore_rng <- sweetnet_restore_rng()
  on.exit(restore_rng(), add = TRUE)
  glycangt_dependencies()
  current <- glycangt_metadata()
  if (!is.null(metadata) && !identical(metadata, current)) {
    cli::cli_abort(
      "GlycanGT model, identifier policy, preprocessing or runtime metadata differs; restore the fitted environment.",
      class = "glyfp_error_backend"
    )
  }
  inputs <- glycangt_inputs(glycans, glycangt_vocabulary())
  weights <- glycangt_load_weights(checkpoint)
  embeddings <- if (fit) {
    list()
  } else {
    glycangt_run(inputs, weights, current$config, current$policy)
  }
  list(checkpoint = checkpoint, metadata = current, embeddings = embeddings)
}

fit_glycangt <- function(glycans, checkpoint, unknown, na_action) {
  if (!identical(unknown, "error")) {
    cli::cli_abort(
      "GlycanGT supports only {.code unknown = 'error'}.",
      class = "glyfp_error_argument"
    )
  }
  response <- glycangt_backend(
    unique(glycans[!is.na(glycans)]),
    checkpoint,
    fit = TRUE
  )
  structure(
    list(
      features = glycangt_features(),
      method = "glycangt",
      params = list(
        unknown = unknown,
        na_action = na_action,
        checkpoint = response$checkpoint,
        model = response$metadata
      ),
      schema_version = 1L,
      versions = engine_versions()
    ),
    class = "glyfp_encoder"
  )
}

transform_glycangt <- function(glycans, encoder) {
  valid <- which(!is.na(glycans))
  unique_glycans <- unique(glycans[valid])
  response <- glycangt_backend(
    unique_glycans,
    encoder$params$checkpoint,
    encoder$params$model
  )
  values <- response$embeddings
  if (
    length(values) != length(unique_glycans) ||
      any(vapply(
        values,
        function(x) length(x) != 768L || any(!is.finite(x)),
        TRUE
      ))
  ) {
    cli::cli_abort(
      "GlycanGT returned invalid embeddings.",
      class = "glyfp_error_backend"
    )
  }
  result <- lapply(seq_along(glycans), function(i) {
    stats::setNames(rep(NA_real_, 768L), encoder$features)
  })
  result[valid] <- lapply(match(glycans[valid], unique_glycans), function(i) {
    stats::setNames(values[[i]], encoder$features)
  })
  names(result) <- names(glycans)
  result
}

valid_glycangt_params <- function(x) {
  is.list(x) &&
    identical(names(x), c("unknown", "na_action", "checkpoint", "model")) &&
    identical(x$unknown, "error") &&
    is.character(x$na_action) &&
    length(x$na_action) == 1L &&
    !is.na(x$na_action) &&
    x$na_action %in% c("propagate", "error") &&
    is.character(x$checkpoint) &&
    length(x$checkpoint) == 1L &&
    !is.na(x$checkpoint) &&
    nzchar(x$checkpoint) &&
    is.list(x$model) &&
    identical(x$model$model, "glycangt-large-mask35") &&
    identical(x$model$config, glycangt_spec()$config) &&
    is.list(x$model$policy)
}
