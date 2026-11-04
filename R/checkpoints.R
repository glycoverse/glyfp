sweetnet_checkpoints <- function() {
  list(
    `sweetnet-v1` = list(
      revision = "0a3981df888532e52b4c2d910bfb56b5a45ad955",
      file = "glycowork_sweetnet_species.pt",
      sha256 = "4ee19fc0593db0562c16585b556050694144816754750ebfe1876ed21ba0f4be"
    ),
    `sweetnet-v2` = list(
      revision = "52187c054b9c988ddb2a5a5af9ee64b51004fd1e",
      file = "glycowork_sweetnet_species_v2.pt",
      sha256 = "bc3d901e4284554c1a245b8e6b7e291a5631c8932ca9c27eb8668b6d4424a146"
    )
  )
}

sweetnet_cache_dir <- function(path) {
  if (
    !is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)
  ) {
    cli::cli_abort(
      "{.arg cache_dir} must be one directory path.",
      class = "glyfp_error_argument"
    )
  }
  path <- path.expand(path)
  if (!grepl("^(/|[[:alpha:]]:|\\\\)", path)) {
    path <- file.path(getwd(), path)
  }
  normalizePath(path, mustWork = FALSE)
}

sweetnet_create_cache <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

sweetnet_download <- function(url, destfile) {
  old <- options(timeout = max(120, getOption("timeout", 60)))
  on.exit(options(old), add = TRUE)
  utils::download.file(
    url,
    destfile,
    method = "libcurl",
    mode = "wb",
    quiet = TRUE
  )
}

sweetnet_resolve <- function(checkpoint, cache_dir, download) {
  checkpoint_resolve(
    checkpoint,
    cache_dir,
    download,
    sweetnet_checkpoints(),
    "SweetNet"
  )
}

checkpoint_resolve <- function(
  checkpoint,
  cache_dir,
  download,
  catalog,
  model
) {
  if (
    !is.character(checkpoint) ||
      length(checkpoint) != 1L ||
      is.na(checkpoint) ||
      !nzchar(checkpoint) ||
      grepl("^[A-Za-z][A-Za-z0-9+.-]*://", checkpoint)
  ) {
    cli::cli_abort(
      "{.arg checkpoint} must be a supported model name ({.or {names(catalog)}}) or a local weight-file path.",
      class = "glyfp_error_argument"
    )
  }
  if (!is.logical(download) || length(download) != 1L || is.na(download)) {
    cli::cli_abort(
      "{.arg download} must be TRUE or FALSE.",
      class = "glyfp_error_argument"
    )
  }
  if (!checkpoint %in% names(catalog)) {
    path <- sweetnet_path(checkpoint, "checkpoint")
    hash <- sweetnet_sha256(path)
    id <- names(catalog)[vapply(
      catalog,
      function(x) identical(x$sha256, hash),
      TRUE
    )]
    if (length(id) != 1L) {
      cli::cli_abort(
        "Local checkpoint SHA-256 does not match a supported official {model} artifact.",
        class = "glyfp_error_checksum"
      )
    }
    return(list(path = path, id = id, source = path))
  }
  id <- checkpoint
  spec <- catalog[[id]]
  cache_dir <- sweetnet_cache_dir(cache_dir)
  path <- file.path(cache_dir, paste0(id, "-", spec$sha256, ".pt"))
  verify <- function(path) {
    file.exists(path) &&
      !dir.exists(path) &&
      identical(sweetnet_sha256(path), spec$sha256)
  }
  if (file.exists(path)) {
    if (!verify(path)) {
      cli::cli_abort(
        "Cached checkpoint failed SHA-256 validation: {.file {path}}. Remove this file to download a fresh copy, or provide a verified local checkpoint.",
        class = "glyfp_error_checksum"
      )
    }
    return(list(path = normalizePath(path), id = id, source = id))
  }
  if (!download) {
    cli::cli_abort(
      "Checkpoint {.val {id}} is not cached at {.file {path}} and {.code download = FALSE}. Supply a local checkpoint or allow its first download.",
      class = "glyfp_error_offline"
    )
  }
  if (
    !dir.exists(cache_dir) &&
      !sweetnet_create_cache(cache_dir) &&
      !dir.exists(cache_dir)
  ) {
    cli::cli_abort(
      "Cannot create checkpoint cache {.file {cache_dir}}.",
      class = "glyfp_error_cache"
    )
  }
  temporary <- tempfile(paste0(id, "-"), tmpdir = cache_dir, fileext = ".part")
  on.exit(unlink(temporary), add = TRUE)
  url <- paste0(
    "https://huggingface.co/DBojar/glycowork_models/resolve/",
    spec$revision,
    "/",
    spec$file
  )
  cli::cli_inform("Downloading {.val {id}} weights to {.file {cache_dir}}.")
  status <- tryCatch(
    suppressWarnings(sweetnet_download(url, temporary)),
    error = function(e) {
      cli::cli_abort(
        "Could not download {.val {id}}. Retry when online or provide a verified local checkpoint.",
        parent = e,
        class = "glyfp_error_download"
      )
    }
  )
  if (!identical(status, 0L)) {
    cli::cli_abort(
      "Download of {.val {id}} did not complete. Retry when online or provide a local checkpoint.",
      class = "glyfp_error_download"
    )
  }
  if (!verify(temporary)) {
    cli::cli_abort(
      "Downloaded checkpoint failed SHA-256 validation; it was not cached.",
      class = "glyfp_error_checksum"
    )
  }
  # Publish only a fully verified file, on the same filesystem. Another process
  # may have populated the cache while this process was downloading.
  if (file.exists(path)) {
    if (!verify(path)) {
      cli::cli_abort(
        "Another process created an invalid cache file at {.file {path}}.",
        class = "glyfp_error_checksum"
      )
    }
  } else if (!file.rename(temporary, path)) {
    if (!verify(path)) {
      cli::cli_abort(
        "Cannot save verified checkpoint at {.file {path}}.",
        class = "glyfp_error_cache"
      )
    }
  }
  list(path = normalizePath(path), id = id, source = id)
}
