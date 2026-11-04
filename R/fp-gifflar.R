#' Encode Glycans with GIFFLAR
#'
#' @description
#' Returns 384-dimensional molecular graph embeddings from the official
#' GIFFLAR species model with verified checkpoint downloads and caching.
#' @inheritParams fp_motif_count
#' @param checkpoint `"gifflar-species"` (default), or the local path to the
#'   official `glycowork_gifflar_species.pt` artifact. The named model downloads
#'   from a pinned official revision on first use and is SHA-256 checked before
#'   caching and loading. Local files are identified by SHA-256, regardless of
#'   filename, and are never replaced or downloaded.
#' @inheritParams fp_sweetnet
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
fp_gifflar <- function(
  glycans,
  checkpoint = "gifflar-species",
  ...,
  cache_dir = tools::R_user_dir("glyfp", "cache"),
  download = TRUE,
  na_action = "propagate"
) {
  rlang::check_dots_empty()
  glycans <- complete_glycans(glycans, na_action)
  embedding_fingerprints(
    glycans,
    run = \(x) {
      gifflar_backend(x, checkpoint, cache_dir = cache_dir, download = download)
    },
    features = gifflar_features(),
    model = "gifflar"
  )
}
