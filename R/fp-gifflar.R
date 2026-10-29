#' Encode Glycans with GIFFLAR
#'
#' @description
#' Returns 384-dimensional molecular graph embeddings from the official
#' GIFFLAR species model using a manually supplied local checkpoint.
#' @inheritParams fp_motif_count
#' @param checkpoint Local path to the official `glycowork_gifflar_species.pt`
#'   artifact. Obtain it from
#'   [DBojar/glycowork_models](https://huggingface.co/DBojar/glycowork_models).
#'   The file is identified by SHA-256, regardless of filename, and is never
#'   replaced or downloaded. A local checkpoint is required on every call.
#'
#' @details
#' GIFFLAR concatenates mean atom, bond and residue representations, using 27
#' fixed residue templates. It rejects unsupported templates,
#' anomeric-to-anomeric bonds and attachments that would move or extend a
#' modification. See the Get started article for supported labels.
#' @inheritSection fp_sweetnet Embedding contract
#' @inheritSection fp_motif_count Complete structures
#'
#' @returns A list of named double vectors with 384 coordinates (`gifflar_001`
#'   to `gifflar_384`), preserving input length, order and names. Missing inputs
#'   return all-NA vectors with the full coordinate set. No model object is returned.
#' @seealso [fp_motif_count()], [fp_sweetnet()], [fp_glycangt()]
#' @examples
#' checkpoint <- Sys.getenv("GLYFP_TEST_GIFFLAR_CHECKPOINT")
#' if (file.exists(checkpoint) && requireNamespace("torch", quietly = TRUE) &&
#'     torch::torch_is_installed()) {
#'   fp <- fp_gifflar("Gal6S(b1-4)GlcNAc(b1-", checkpoint = checkpoint)
#'   lengths(fp)
#' }
#' @export
fp_gifflar <- function(glycans, checkpoint, ..., na_action = "propagate") {
  rlang::check_dots_empty()
  if (missing(checkpoint)) {
    cli::cli_abort(
      "GIFFLAR requires an explicit local {.arg checkpoint}; download the official weights yourself. glyfp never downloads GIFFLAR weights.",
      class = "glyfp_error_argument"
    )
  }
  glycans <- complete_glycans(glycans, na_action)
  embedding_fingerprints(
    glycans,
    run = \(x) gifflar_backend(x, checkpoint),
    features = gifflar_features(),
    model = "gifflar"
  )
}
