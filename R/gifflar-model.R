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

gifflar_checkpoints <- function() {
  spec <- gifflar_spec()
  list(
    `gifflar-species` = list(
      revision = spec$checkpoint_revision,
      file = "glycowork_gifflar_species.pt",
      sha256 = spec$sha256
    )
  )
}

gifflar_resolve <- function(checkpoint, cache_dir, download) {
  checkpoint_resolve(
    checkpoint,
    cache_dir,
    download,
    gifflar_checkpoints(),
    "GIFFLAR"
  )
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

gifflar_backend <- function(
  glycans,
  checkpoint,
  cache_dir = tools::R_user_dir("glyfp", "cache"),
  download = TRUE
) {
  checkpoint <- gifflar_resolve(checkpoint, cache_dir, download)$path
  restore_rng <- sweetnet_restore_rng()
  on.exit(restore_rng(), add = TRUE)
  gifflar_dependencies()
  inputs <- gifflar_inputs(glycans, gifflar_templates(), sweetnet_vocabulary())
  weights <- gifflar_load_weights(checkpoint)
  lapply(inputs, function(x) gifflar_infer(x, weights))
}
