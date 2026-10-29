# Native GIFFLAR network, adapted from MIT-licensed glycowork commit
# 199efdde6bc137dbbcd900f84f927fdbbca52eae. See inst/COPYRIGHTS.
gifflar_relations <- function() {
  c(
    "atoms___coboundary___atoms",
    "atoms___to___bonds",
    "bonds___to___monosacchs",
    "bonds___boundary___bonds",
    "bonds___coboundary___bonds",
    "monosacchs___boundary___monosacchs"
  )
}

gifflar_infer <- function(input, weights, layers = 4L) {
  torch::with_no_grad({
    index <- function(x) {
      torch::torch_tensor(
        as.integer(x) + 1L,
        dtype = torch::torch_long(),
        device = "cpu"
      )
    }
    embedding <- function(ids, key) weights[[key]]$index_select(1L, index(ids))
    linear <- function(x, prefix) {
      torch::nnf_linear(
        x,
        weights[[paste0(prefix, ".weight")]],
        weights[[paste0(prefix, ".bias")]]
      )
    }
    norm <- function(x, prefix) {
      torch::nnf_batch_norm(
        x,
        running_mean = weights[[paste0(prefix, ".running_mean")]],
        running_var = weights[[paste0(prefix, ".running_var")]],
        weight = weights[[paste0(prefix, ".weight")]],
        bias = weights[[paste0(prefix, ".bias")]],
        training = FALSE,
        eps = 1e-5
      )
    }
    activate <- function(x, prefix) {
      torch::nnf_prelu(x, weights[[paste0(prefix, ".weight")]])
    }
    types <- c("atoms", "bonds", "monosacchs")
    keys <- c(
      "atom_embedding.weight",
      "bond_embedding.weight",
      "mono_embedding.weight"
    )
    x <- stats::setNames(
      lapply(seq_along(types), function(i) {
        embedding(input$x[[types[i]]], keys[i])
      }),
      types
    )
    for (layer in seq_len(layers) - 1L) {
      outputs <- stats::setNames(vector("list", length(types)), types)
      for (relation in gifflar_relations()) {
        edges <- input$edges[[relation]]
        if (is.null(edges) || !nrow(edges)) {
          next
        }
        parts <- strsplit(relation, "___", fixed = TRUE)[[1L]]
        src <- x[[parts[1L]]]
        dst <- x[[parts[3L]]]
        prefix <- paste0("convs.", layer, ".convs.<", relation, ">")
        messages <- src$index_select(1L, index(edges[, 1L]))
        if (identical(relation, "monosacchs___boundary___monosacchs")) {
          attr <- embedding(
            input$edge_attr[[relation]],
            "link_embedding.weight"
          )
          messages <- torch::nnf_relu(
            messages + linear(attr, paste0(prefix, ".lin"))
          )
        }
        summed <- torch::torch_zeros_like(dst)
        summed$index_add_(1L, index(edges[, 2L]), messages)
        y <- summed + (1 + weights[[paste0(prefix, ".eps")]]) * dst
        y <- linear(
          activate(
            norm(linear(y, paste0(prefix, ".nn.0")), paste0(prefix, ".nn.1")),
            paste0(prefix, ".nn.2")
          ),
          paste0(prefix, ".nn.3")
        )
        outputs[[parts[3L]]] <- c(outputs[[parts[3L]]], list(y))
      }
      for (type in types) {
        values <- outputs[[type]]
        if (!length(values)) {
          cli::cli_abort(
            "GIFFLAR graph lacks incoming relations for {.val {type}}.",
            class = "glyfp_error_backend"
          )
        }
        y <- torch::torch_stack(values, dim = 1L)$sum(dim = 1L)
        y <- activate(
          norm(y, paste0("norms.", layer, ".", type)),
          paste0("acts.", layer, ".", type)
        )
        x[[type]] <- y + x[[type]]
      }
    }
    as.double(torch::as_array(torch::torch_cat(
      lapply(x, function(y) y$mean(dim = 1L)),
      dim = 1L
    )))
  })
}
