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

glycangt_backend <- function(glycans, checkpoint) {
  checkpoint <- glycangt_checkpoint(checkpoint)
  restore_rng <- sweetnet_restore_rng()
  on.exit(restore_rng(), add = TRUE)
  glycangt_dependencies()
  inputs <- glycangt_inputs(glycans, glycangt_vocabulary())
  weights <- glycangt_load_weights(checkpoint)
  glycangt_run(inputs, weights, glycangt_spec()$config, glycangt_policy())
}
