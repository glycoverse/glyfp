# Version 1 uses byte-length-prefixed UTF-8 fields, never R serialization.
residue_record <- function(fields) {
  fields <- enc2utf8(as.character(fields))
  paste0(nchar(fields, type = "bytes"), ":", fields, collapse = "")
}

residue_hash <- function(fields) {
  digest::digest(residue_record(fields), algo = "sha256", serialize = FALSE)
}

validate_residue_options <- function(radius, n_bits, mode) {
  integer_in <- function(x, lower, upper) {
    is.numeric(x) &&
      length(x) == 1L &&
      !is.na(x) &&
      is.finite(x) &&
      x >= lower &&
      x <= upper &&
      x == floor(x)
  }
  if (!integer_in(radius, 0, 10)) {
    cli::cli_abort(
      "{.arg radius} must be an integer from 0 to 10.",
      class = "glyfp_error_argument"
    )
  }
  if (!integer_in(n_bits, 1, 1048576)) {
    cli::cli_abort(
      "{.arg n_bits} must be an integer from 1 to 1048576.",
      class = "glyfp_error_argument"
    )
  }
  validate_choice(mode, c("binary", "count"), "mode")
}

# A node's identity includes represented residue stereochemistry, substituents,
# its own anomer/donor carbon and whether it is the reducing-end root. Edge
# messages include the full bond label and donor/acceptor direction.
residue_environments <- function(graph, radius) {
  n <- igraph::vcount(graph)
  edges <- igraph::as_edgelist(graph, names = FALSE)
  bonds <- igraph::edge_attr(graph, "linkage")
  root <- which(igraph::degree(graph, mode = "in") == 0L)
  anomer <- character(n)
  anomer[root] <- igraph::graph_attr(graph, "anomer")
  if (nrow(edges)) {
    anomer[edges[, 2L]] <- sub("-.*$", "", bonds)
  }
  mono <- igraph::vertex_attr(graph, "mono")
  subs <- igraph::vertex_attr(graph, "sub")
  labels <- vapply(
    seq_len(n),
    \(i) {
      residue_hash(c(
        "glyfp-residue-v1",
        mono[i],
        subs[i],
        anomer[i],
        if (i %in% root) "root" else "internal"
      ))
    },
    ""
  )
  rounds <- list(labels)
  neighbors <- vector("list", n)
  for (e in seq_len(nrow(edges))) {
    a <- edges[e, 1L]
    b <- edges[e, 2L]
    neighbors[[a]] <- c(
      neighbors[[a]],
      list(list(node = b, role = "acceptor", bond = bonds[e]))
    )
    neighbors[[b]] <- c(
      neighbors[[b]],
      list(list(node = a, role = "donor", bond = bonds[e]))
    )
  }
  for (iteration in seq_len(radius)) {
    labels <- vapply(
      seq_len(n),
      \(i) {
        messages <- vapply(
          neighbors[[i]],
          \(edge) {
            residue_record(c(edge$role, edge$bond, labels[edge$node]))
          },
          ""
        )
        residue_hash(c(
          "glyfp-wl-v1",
          iteration,
          labels[i],
          sort(messages, method = "radix")
        ))
      },
      ""
    )
    rounds[[iteration + 1L]] <- labels
  }
  rounds
}

residue_tokens <- function(graph, radius, method) {
  rounds <- residue_environments(graph, radius)
  if (method == "circular") {
    return(unlist(
      lapply(seq_along(rounds), \(k) {
        vapply(
          rounds[[k]],
          \(label) residue_hash(c("glyfp-circular-v1", k - 1L, label)),
          "",
          USE.NAMES = FALSE
        )
      }),
      use.names = FALSE
    ))
  }
  n <- igraph::vcount(graph)
  if (n < 2L) {
    return(character())
  }
  labels <- rounds[[radius + 1L]]
  distances <- igraph::distances(graph, mode = "all", weights = NA)
  # Every distinct unordered pair contributes once. There are no self-pairs.
  unlist(
    lapply(seq_len(n - 1L), \(i) {
      vapply(
        seq.int(i + 1L, n),
        \(j) {
          residue_hash(c(
            "glyfp-residue-pair-v1",
            as.integer(distances[i, j]),
            sort(c(labels[i], labels[j]), method = "radix")
          ))
        },
        ""
      )
    }),
    use.names = FALSE
  )
}

residue_bins <- function(tokens, n_bits) {
  # Big-endian first 32 SHA-256 bits modulo n_bits. Horner reduction avoids
  # signed-integer overflow and is exact in R doubles for all allowed n_bits.
  vapply(
    substr(tokens, 1L, 8L),
    \(token) {
      digits <- strsplit(token, "", fixed = TRUE)[[1L]]
      value <- 0
      for (digit in digits) {
        value <- (16 * value + strtoi(digit, 16L)) %% n_bits
      }
      as.integer(value + 1)
    },
    1L,
    USE.NAMES = FALSE
  )
}

hashed_residue_fingerprints <- function(glycans, radius, n_bits, mode, method) {
  features <- sprintf(
    "%s_v1_r%d_n%d_%s_%07d",
    method,
    radius,
    n_bits,
    mode,
    seq_len(n_bits)
  )
  valid <- which(!is.na(glycans))
  unique_glycans <- unique(glycans[valid])
  values <- lapply(as.list(unique_glycans), \(graph) {
    tokens <- residue_tokens(graph, radius, method)
    counts <- tabulate(residue_bins(tokens, n_bits), nbins = n_bits)
    if (mode == "binary") {
      counts <- as.integer(counts > 0L)
    }
    stats::setNames(counts, features)
  })
  result <- lapply(seq_along(glycans), \(i) {
    stats::setNames(rep(NA_integer_, n_bits), features)
  })
  result[valid] <- values[match(glycans[valid], unique_glycans)]
  names(result) <- names(glycans)
  result
}
