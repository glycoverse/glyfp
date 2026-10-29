# Internal GlycanGT inference building blocks. Not registered as a public method
# until the gated pretrained checkpoint has passed numerical validation.
# Adapted under Apache-2.0; see inst/licenses/GlycanGT-LICENSE.
# Architecture: matsui-lab/GlycanGT 96611518c971deb89215ca163deaf9de3a59fa32.

glycangt_policy <- function(seed = 1729L, dimensions = 64L) {
  valid_integer <- function(x) {
    is.numeric(x) &&
      length(x) == 1L &&
      !is.na(x) &&
      is.finite(x) &&
      x >= 0 &&
      x <= .Machine$integer.max &&
      x == as.integer(x)
  }
  if (
    !valid_integer(seed) ||
      seed < 0 ||
      seed > .Machine$integer.max ||
      !valid_integer(dimensions) ||
      dimensions < 1L
  ) {
    cli::cli_abort(
      "Invalid GlycanGT random identifier policy.",
      class = "glyfp_error_argument"
    )
  }
  list(
    version = 1L,
    algorithm = "sha256-canonical-seed7hex:torch-cpu-gaussian-qr",
    seed = as.integer(seed),
    dimensions = as.integer(dimensions),
    dtype = "float32",
    graph_order = "glyrepr-canonical-iupac",
    batch_size = 1L
  )
}

glycangt_validate_policy <- function(policy) {
  if (
    !is.list(policy) ||
      !identical(policy, glycangt_policy(policy$seed, policy$dimensions))
  ) {
    cli::cli_abort(
      "Unsupported GlycanGT random identifier policy.",
      class = "glyfp_error_backend"
    )
  }
  invisible(policy)
}

glycangt_graph_seed <- function(key, policy) {
  glycangt_validate_policy(policy)
  # A stable independent stream per canonical graph: duplicates and batch order
  # do not consume the caller's stream or change the graph's random identifiers.
  hash <- digest::digest(
    paste0(policy$seed, "\n", enc2utf8(key)),
    algo = "sha256",
    serialize = FALSE
  )
  strtoi(substr(hash, 1L, 7L), base = 16L)
}

glycangt_orf <- function(nodes, key, policy) {
  restore_rng <- sweetnet_restore_rng()
  on.exit(restore_rng(), add = TRUE)
  seed <- glycangt_graph_seed(key, policy)
  # This is the official one-graph distribution: Gaussian n x n -> QR transpose
  # -> row normalization -> pad/truncate to Do -> row normalization. Do not
  # replace it with a Do x Do basis or a fixed identity matrix.
  torch::with_torch_manual_seed(
    {
      gaussian <- torch::torch_randn(
        c(1L, nodes, nodes),
        device = "cpu",
        dtype = torch::torch_float32()
      )
      q <- torch::linalg_qr(gaussian, mode = "reduced")[[1L]]$transpose(2L, 3L)
      q <- torch::nnf_normalize(q, p = 2, dim = 3L)$squeeze(1L)
      if (nodes < policy$dimensions) {
        q <- torch::nnf_pad(q, c(0L, policy$dimensions - nodes))
      } else if (nodes > policy$dimensions) {
        q <- q[, seq_len(policy$dimensions), drop = FALSE]
      }
      torch::nnf_normalize(q, p = 2, dim = 2L)
    },
    seed = seed
  )
}

glycangt_inputs <- function(glycans, vocabulary) {
  glycans <- normalize_glycans(glycans, "propagate")
  result <- lapply(as.list(glycans), function(g) {
    if (is.null(g)) {
      return(NULL)
    }
    # Reconstruct through the public constructor to canonicalize graph storage.
    structure <- glyrepr::glycan_structure(g)
    g <- as.list(structure)[[1L]]
    mono <- igraph::vertex_attr(g, "mono")
    root <- which(igraph::degree(g, mode = "in") == 0L)
    expected <- paste0("?", glyrepr::infer_anomer_pos(mono[root]))
    if (
      length(mono) < 2L ||
        length(root) != 1L ||
        !identical(igraph::graph_attr(g, "anomer"), expected)
    ) {
      cli::cli_abort(
        "GlycanGT currently requires at least two residues and an unspecified normal reducing end.",
        class = "glyfp_error_unsupported"
      )
    }
    linkage <- igraph::edge_attr(g, "linkage")
    if (any(!grepl("^[ab][0-9]+-[0-9]+$", linkage))) {
      cli::cli_abort(
        "GlycanGT currently requires fully specified internal linkages.",
        class = "glyfp_error_unsupported"
      )
    }
    sub <- igraph::vertex_attr(g, "sub")
    sub[is.na(sub)] <- ""
    labels <- paste0(mono, gsub(",", "", sub, fixed = TRUE))
    node_ids <- unname(vocabulary$monomer[labels])
    edge_ids <- unname(vocabulary$linkage[linkage])
    if (anyNA(c(node_ids, edge_ids))) {
      cli::cli_abort(
        "GlycanGT input contains a token outside the supplied vocabulary.",
        class = "glyfp_error_unknown"
      )
    }
    # Official triples are parent -> child; GlycanGraphDataset reverses them.
    endpoints <- igraph::as_edgelist(g, names = FALSE)
    list(
      key = unname(as.character(structure)),
      nodes = as.integer(node_ids),
      edges = as.integer(edge_ids),
      index = unname(endpoints[, 2:1, drop = FALSE] - 1L)
    )
  })
  names(result) <- names(glycans)
  result
}

glycangt_infer <- function(input, weights, config, policy) {
  glycangt_validate_policy(policy)
  if (
    policy$dimensions != config$orf_dimensions ||
      config$dimensions %% config$heads != 0L
  ) {
    cli::cli_abort(
      "GlycanGT architecture and identifier policy differ.",
      class = "glyfp_error_backend"
    )
  }
  torch::with_no_grad({
    dimensions <- config$dimensions
    n <- length(input$nodes)
    e <- length(input$edges)
    index <- function(x) {
      torch::torch_tensor(x + 1L, dtype = torch::torch_long(), device = "cpu")
    }
    embedding <- function(ids, name) {
      torch::nnf_embedding(index(ids), weights[[name]])
    }
    prefix <- "encoder.graph_feature."
    x <- torch::torch_cat(
      list(
        embedding(input$nodes, paste0(prefix, "atom_encoder.weight")),
        embedding(input$edges, paste0(prefix, "edge_encoder.weight"))
      ),
      dim = 1L
    )
    identifiers <- glycangt_orf(n, input$key, policy)
    pairs <- rbind(cbind(seq_len(n) - 1L, seq_len(n) - 1L), input$index)
    pair_features <- torch::torch_cat(
      list(
        identifiers$index_select(1L, index(pairs[, 1L])),
        identifiers$index_select(1L, index(pairs[, 2L]))
      ),
      dim = 2L
    )
    x <- x +
      torch::nnf_linear(
        pair_features,
        weights[[paste0(prefix, "orf_encoder.weight")]]
      ) +
      embedding(
        as.integer(pairs[, 1L] == pairs[, 2L]),
        paste0(prefix, "order_encoder.weight")
      )
    x <- torch::torch_cat(
      list(
        weights[[paste0(prefix, "graph_token.weight")]],
        weights[[paste0(prefix, "null_token.weight")]],
        x
      ),
      dim = 1L
    )
    tokens <- n + e + 2L
    head_size <- dimensions %/% config$heads
    for (layer in seq_len(config$layers)) {
      base <- paste0("encoder.layers.", layer - 1L, ".")
      linear <- function(y, name) {
        torch::nnf_linear(
          y,
          weights[[paste0(base, name, ".weight")]],
          weights[[paste0(base, name, ".bias")]]
        )
      }
      norm <- function(y, name) {
        torch::nnf_layer_norm(
          y,
          dimensions,
          weights[[paste0(base, name, ".weight")]],
          weights[[paste0(base, name, ".bias")]],
          eps = 1e-5
        )
      }
      split_heads <- function(y) {
        y$reshape(c(tokens, config$heads, head_size))$transpose(1L, 2L)
      }
      q <- split_heads(linear(x, "self_attn.q_proj") * head_size^(-0.5))
      k <- split_heads(linear(x, "self_attn.k_proj"))
      v <- split_heads(linear(x, "self_attn.v_proj"))
      attention <- torch::nnf_softmax(
        torch::torch_bmm(q, k$transpose(2L, 3L)),
        dim = 3L
      )
      y <- torch::torch_bmm(attention, v)$transpose(
        1L,
        2L
      )$contiguous()$reshape(c(tokens, dimensions))
      x <- norm(x + linear(y, "self_attn.out_proj"), "self_attn_layer_norm")
      y <- linear(
        torch::nnf_gelu(linear(x, "feedforward.fc1"), approximate = "none"),
        "feedforward.fc2"
      )
      x <- norm(x + y, "final_layer_norm")
    }
    as.double(torch::as_array(x[1L, ]))
  })
}

glycangt_run <- function(inputs, weights, config, policy) {
  restore_rng <- sweetnet_restore_rng()
  on.exit(restore_rng(), add = TRUE)
  # Weights are explicitly supplied by internal callers; this function neither
  # downloads nor constructs a model. The pretrained adapter is pending access.
  glycangt_validate_policy(policy)
  glycangt_validate_weights(weights, config)
  lapply(inputs, function(input) {
    if (is.null(input)) {
      return(rep(NA_real_, config$dimensions))
    }
    glycangt_infer(input, weights, config, policy)
  })
}


glycangt_validate_weights <- function(weights, config) {
  fields <- c("dimensions", "heads", "layers", "orf_dimensions")
  invalid <- !is.list(config) ||
    !identical(names(config), fields) ||
    !all(vapply(
      config,
      function(x) {
        is.numeric(x) &&
          length(x) == 1L &&
          is.finite(x) &&
          x >= 1 &&
          x <= .Machine$integer.max &&
          x == as.integer(x)
      },
      TRUE
    ))
  if (invalid || config$dimensions %% config$heads != 0L || !is.list(weights)) {
    cli::cli_abort(
      "Unsupported GlycanGT inference configuration.",
      class = "glyfp_error_backend"
    )
  }
  d <- config$dimensions
  check <- function(name, shape) {
    x <- weights[[name]]
    if (
      !inherits(x, "torch_tensor") ||
        !identical(as.integer(x$size()), as.integer(shape))
    ) {
      cli::cli_abort(
        "Unexpected GlycanGT tensor shape for {.val {name}}.",
        class = "glyfp_error_backend"
      )
    }
  }
  prefix <- "encoder.graph_feature."
  for (name in c("atom_encoder", "edge_encoder")) {
    w <- weights[[paste0(prefix, name, ".weight")]]
    if (
      !inherits(w, "torch_tensor") || length(w$size()) != 2L || w$size(1L) < 2L
    ) {
      cli::cli_abort(
        "Invalid GlycanGT vocabulary embedding.",
        class = "glyfp_error_backend"
      )
    }
    check(paste0(prefix, name, ".weight"), c(w$size(1L), d))
  }
  for (name in c("graph_token", "null_token")) {
    check(paste0(prefix, name, ".weight"), c(1L, d))
  }
  check(paste0(prefix, "orf_encoder.weight"), c(d, 2L * config$orf_dimensions))
  check(paste0(prefix, "order_encoder.weight"), c(2L, d))
  for (layer in seq_len(config$layers)) {
    base <- paste0("encoder.layers.", layer - 1L, ".")
    for (name in c("q_proj", "k_proj", "v_proj", "out_proj")) {
      check(paste0(base, "self_attn.", name, ".weight"), c(d, d))
      check(paste0(base, "self_attn.", name, ".bias"), d)
    }
    for (name in c("self_attn_layer_norm", "final_layer_norm")) {
      check(paste0(base, name, ".weight"), d)
      check(paste0(base, name, ".bias"), d)
    }
    w <- weights[[paste0(base, "feedforward.fc1.weight")]]
    if (!inherits(w, "torch_tensor") || length(w$size()) != 2L) {
      cli::cli_abort(
        "Invalid GlycanGT feedforward tensor.",
        class = "glyfp_error_backend"
      )
    }
    f <- w$size(1L)
    check(paste0(base, "feedforward.fc1.weight"), c(f, d))
    check(paste0(base, "feedforward.fc1.bias"), f)
    check(paste0(base, "feedforward.fc2.weight"), c(d, f))
    check(paste0(base, "feedforward.fc2.bias"), d)
  }
  invisible(TRUE)
}
