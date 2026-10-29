#' Encode Glycans with GlycanGT
#'
#' @description
#' Returns 768-dimensional graph Transformer embeddings from the official
#' GlycanGT large mask-35 model using a manually supplied local checkpoint.
#' @inheritParams fp_motif_count
#' @param checkpoint Local path to the official
#'   `pretrained_tokengt_large_multitask_mask_35_final.pt` artifact. Obtain it from
#'   [Akikitani295/GlycanGT](https://huggingface.co/Akikitani295/GlycanGT).
#'   The file is identified by SHA-256, regardless of filename, and is never
#'   replaced or downloaded. A local checkpoint is required on every call.
#'
#' @details
#' GlycanGT returns the model's final graph token and requires at least two
#' residues. Its orthogonal random identifiers use master seed 1729, a
#' SHA-256-derived seed per canonical projected glycan, the official single-graph
#' Gaussian QR construction and width 64. This fixed policy preserves caller R
#' and CPU torch RNG state. Upstream multi-graph ORF batches can yield different
#' outputs; this adapter does not claim to reproduce published downstream performance.
#' @inheritSection fp_sweetnet Embedding contract
#' @inheritSection fp_motif_count Complete structures
#'
#' @returns A list of named double vectors with 768 coordinates (`glycangt_001`
#'   to `glycangt_768`), preserving input length, order and names. Missing inputs
#'   return all-NA vectors with the full coordinate set. No model object is returned.
#' @seealso [fp_motif_count()], [fp_sweetnet()], [fp_gifflar()]
#' @examples
#' checkpoint <- Sys.getenv("GLYFP_TEST_GLYCANGT_CHECKPOINT")
#' if (file.exists(checkpoint) && requireNamespace("torch", quietly = TRUE) &&
#'     torch::torch_is_installed()) {
#'   fp <- fp_glycangt("Gal(b1-4)GlcNAc(b1-", checkpoint = checkpoint)
#'   lengths(fp)
#' }
#' @export
fp_glycangt <- function(glycans, checkpoint, ..., na_action = "propagate") {
  rlang::check_dots_empty()
  if (missing(checkpoint)) {
    cli::cli_abort(
      "GlycanGT requires an explicit local {.arg checkpoint}; download the official weights yourself. glyfp never downloads GlycanGT weights.",
      class = "glyfp_error_argument"
    )
  }
  glycans <- complete_glycans(glycans, na_action)
  embedding_fingerprints(
    glycans,
    run = \(x) glycangt_backend(x, checkpoint),
    features = glycangt_features(),
    model = "glycangt"
  )
}
