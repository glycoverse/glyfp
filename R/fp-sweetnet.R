#' Encode Glycans with SweetNet
#'
#' @description
#' Returns pretrained SweetNet embeddings, with a choice of official model
#' versions and control over checkpoint caching and downloads.
#' @inheritParams fp_motif_count
#' @param checkpoint `"sweetnet-v2"` (default), `"sweetnet-v1"`, or the local
#'   path to either supported official artifact. Named models download from a
#'   pinned official revision on first use and are SHA-256 checked before caching
#'   and loading. Local files are identified by SHA-256, regardless of filename,
#'   and are never replaced or downloaded.
#' @param cache_dir Named-checkpoint cache directory. No directory is created
#'   for local checkpoints or offline cache misses.
#' @param download Whether to download a missing named checkpoint. Set `FALSE`
#'   for offline use.
#'
#' @details
#' SweetNet v1 and v2 return 128 coordinates from `bn2(lin2(...))`, before final
#' activation, dropout and classification. The default v2 follows the upstream
#' species model; v1 and v2 coordinate spaces must not be mixed.
#'
#' @section Embedding contract:
#' Pretrained models do not encode reducing-end alpha/beta information: only
#' after complete-structure validation, preprocessing omits this known state at
#' the normal reducing carbon. Thus embeddings cannot distinguish glycans
#' differing only in this state; motif fingerprints retain it. Other reducing-end
#' positions and unsupported chemistry are rejected. Tokens outside the fixed
#' model vocabulary always raise an error.
#'
#' Embeddings run in R with optional torch and digest; SweetNet and GlycanGT also
#' use zip. No Python, RDKit or GlyLES runtime is used. Install torch's runtime
#' explicitly: loading glyfp and calling inference never install it automatically.
#' CPU evaluation uses float32, frozen normalization, no dropout or gradients,
#' and one glycan at a time. Batch order, duplicates and batch size do not change
#' coordinates under the same model, weights, preprocessing and runtime settings.
#' Preserve these settings and package/runtime versions with saved results.
#' Cross-platform comparisons should allow floating-point tolerance (the reference
#' tests bound maximum absolute coordinate error by 1e-5).
#' @inheritSection fp_motif_count Complete structures
#'
#' @returns A list of named double vectors with 128 coordinates (`sweetnet_001`
#'   to `sweetnet_128`), preserving input length, order and names. Missing inputs
#'   return all-NA vectors with the full coordinate set. No model object is returned.
#' @seealso [fp_motif_count()], [fp_glycangt()], [fp_gifflar()]
#' @examples
#' # Run when an official checkpoint is already cached; no download is allowed.
#' cache <- tools::R_user_dir("glyfp", "cache")
#' if (requireNamespace("torch", quietly = TRUE) &&
#'     torch::torch_is_installed() &&
#'     length(list.files(cache, "^sweetnet-v2-.*[.]pt$")) > 0L) {
#'   fp <- fp_sweetnet("Gal(b1-4)GlcNAc(b1-", download = FALSE)
#'   lengths(fp)
#' }
#' @export
fp_sweetnet <- function(
  glycans,
  ...,
  checkpoint = "sweetnet-v2",
  cache_dir = tools::R_user_dir("glyfp", "cache"),
  download = TRUE,
  na_action = "propagate"
) {
  rlang::check_dots_empty()
  glycans <- complete_glycans(glycans, na_action)
  embedding_fingerprints(
    glycans,
    run = \(x) {
      sweetnet_run(
        sweetnet_inputs(x),
        checkpoint,
        cache_dir = cache_dir,
        download = download
      )
    },
    features = sweetnet_features(),
    model = "sweetnet"
  )
}
