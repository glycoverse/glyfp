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
  cache_dir = tools::R_user_dir("glyfp", "cache"),
  download = TRUE
) {
  if (!is.logical(download) || length(download) != 1L || is.na(download)) {
    cli::cli_abort(
      "{.arg download} must be TRUE or FALSE.",
      class = "glyfp_error_argument"
    )
  }
  cache_dir <- sweetnet_cache_dir(cache_dir)
  restore_rng <- sweetnet_restore_rng()
  on.exit(restore_rng(), add = TRUE)
  sweetnet_dependencies()
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
  weights <- sweetnet_load_weights(resolved$path, resolved$id)
  lapply(
    inputs,
    sweetnet_infer,
    weights = weights,
    vocabulary = vocabulary,
    id = resolved$id
  )
}
