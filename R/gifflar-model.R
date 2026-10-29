gifflar_spec <- function() {
  list(
    model = "gifflar-species",
    dimensions = 384L,
    layers = 4L,
    source_revision = "199efdde6bc137dbbcd900f84f927fdbbca52eae",
    checkpoint_revision = "52187c054b9c988ddb2a5a5af9ee64b51004fd1e",
    sha256 = "0dab7c3b61aad95dab37d4974da4b57b3fd985128dc09ad5da846e82adb9519d",
    vocabulary_sha256 = "9edf1dbdb4df99ae6444108ba09e72f786c175c0b051a32fd82af3e46edc69fc",
    templates_sha256 = "9e6855f7e2c185076dd65dcf655ebc5b52f68d27b37d8b3d55b6b3bade222acd"
  )
}

gifflar_features <- function() sprintf("gifflar_%03d", seq_len(384L))

gifflar_checkpoint <- function(checkpoint) {
  if (
    !is.character(checkpoint) ||
      length(checkpoint) != 1L ||
      is.na(checkpoint) ||
      !nzchar(checkpoint) ||
      grepl("^[A-Za-z][A-Za-z0-9+.-]*://", checkpoint)
  ) {
    cli::cli_abort(
      "GIFFLAR requires a local {.arg checkpoint} file downloaded from https://huggingface.co/DBojar/glycowork_models. No automatic download is performed.",
      class = "glyfp_error_argument"
    )
  }
  if (!file.exists(checkpoint) || dir.exists(checkpoint)) {
    cli::cli_abort(
      "GIFFLAR checkpoint file does not exist: {.path {checkpoint}}.",
      class = "glyfp_error_backend"
    )
  }
  normalizePath(checkpoint, mustWork = TRUE)
}

gifflar_dependencies <- function() {
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
  for (pkg in c("torch", "digest")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cli::cli_abort(
        "Install {.pkg {pkg}} to use GIFFLAR.",
        class = "glyfp_error_backend"
      )
    }
  }
  if (!torch::torch_is_installed()) {
    cli::cli_abort(
      "Install the R torch runtime with {.code torch::install_torch()} to use GIFFLAR.",
      class = "glyfp_error_backend"
    )
  }
}

.gifflar_state <- new.env(parent = emptyenv())

gifflar_load_weights <- function(checkpoint) {
  if (!identical(sweetnet_sha256(checkpoint), gifflar_spec()$sha256)) {
    cli::cli_abort(
      "Checkpoint SHA-256 does not match the supported official GIFFLAR species artifact.",
      class = "glyfp_error_backend"
    )
  }
  if (is.null(.gifflar_state$weights)) {
    # The pinned official state dictionary is already CPU-tagged. Never alter it.
    .gifflar_state$weights <- torch::load_state_dict(checkpoint)
  }
  .gifflar_state$weights
}

gifflar_metadata <- function() {
  sources <- list(
    spec = gifflar_spec,
    dependencies = gifflar_dependencies,
    rng = sweetnet_restore_rng,
    inputs = gifflar_inputs,
    graph = gifflar_graph,
    smiles = gifflar_smiles,
    parsing = gifflar_parse_smiles,
    templates = gifflar_templates,
    vocabulary = sweetnet_vocabulary,
    relations = gifflar_relations,
    checkpoint = gifflar_checkpoint,
    weights = gifflar_load_weights,
    inference = gifflar_infer
  )
  list(
    model = gifflar_spec()$model,
    spec = gifflar_spec(),
    preprocessing = "glyrepr-canonical:27-fixed-residue-templates:concrete-links:unspecified-root",
    readout = "concat-mean-atoms-bonds-monosacchs",
    dtype = "float32",
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

gifflar_backend <- function(glycans, checkpoint, metadata = NULL, fit = FALSE) {
  checkpoint <- gifflar_checkpoint(checkpoint)
  restore_rng <- sweetnet_restore_rng()
  on.exit(restore_rng(), add = TRUE)
  gifflar_dependencies()
  current <- gifflar_metadata()
  if (!is.null(metadata) && !identical(metadata, current)) {
    cli::cli_abort(
      "GIFFLAR model, preprocessing or runtime metadata differs; restore the fitted environment.",
      class = "glyfp_error_backend"
    )
  }
  inputs <- gifflar_inputs(glycans, gifflar_templates(), sweetnet_vocabulary())
  weights <- gifflar_load_weights(checkpoint)
  embeddings <- if (fit) {
    list()
  } else {
    lapply(inputs, function(x) gifflar_infer(x, weights))
  }
  list(checkpoint = checkpoint, metadata = current, embeddings = embeddings)
}

fit_gifflar <- function(glycans, checkpoint, unknown, na_action) {
  if (!identical(unknown, "error")) {
    cli::cli_abort(
      "GIFFLAR supports only {.code unknown = 'error'}.",
      class = "glyfp_error_argument"
    )
  }
  response <- gifflar_backend(
    unique(glycans[!is.na(glycans)]),
    checkpoint,
    fit = TRUE
  )
  structure(
    list(
      features = gifflar_features(),
      method = "gifflar",
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

transform_gifflar <- function(glycans, encoder) {
  valid <- which(!is.na(glycans))
  unique_glycans <- unique(glycans[valid])
  response <- gifflar_backend(
    unique_glycans,
    encoder$params$checkpoint,
    encoder$params$model
  )
  values <- response$embeddings
  if (
    length(values) != length(unique_glycans) ||
      any(vapply(
        values,
        function(x) length(x) != 384L || any(!is.finite(x)),
        TRUE
      ))
  ) {
    cli::cli_abort(
      "GIFFLAR returned invalid embeddings.",
      class = "glyfp_error_backend"
    )
  }
  result <- lapply(seq_along(glycans), function(i) {
    stats::setNames(rep(NA_real_, 384L), encoder$features)
  })
  result[valid] <- lapply(match(glycans[valid], unique_glycans), function(i) {
    stats::setNames(values[[i]], encoder$features)
  })
  names(result) <- names(glycans)
  result
}

valid_gifflar_params <- function(x) {
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
    identical(x$model$model, "gifflar-species") &&
    identical(x$model$spec, gifflar_spec())
}
