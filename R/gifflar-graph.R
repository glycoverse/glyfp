# Fixed, explicitly supported residue templates preserve the upstream SMILES
# traversal and chirality tags. No chemistry toolkit or Python is used at runtime.
gifflar_templates <- function() {
  path <- system.file(
    "extdata",
    "gifflar-templates.tsv",
    package = "glyfp",
    mustWork = TRUE
  )
  if (!identical(sweetnet_sha256(path), gifflar_spec()$templates_sha256)) {
    cli::cli_abort(
      "GIFFLAR residue templates have changed.",
      class = "glyfp_error_backend"
    )
  }
  utils::read.delim(
    path,
    header = FALSE,
    quote = "",
    comment.char = "",
    colClasses = "character"
  )
}

gifflar_smiles <- function(g, templates) {
  mono <- igraph::vertex_attr(g, "mono")
  sub <- igraph::vertex_attr(g, "sub")
  sub[is.na(sub)] <- ""
  labels <- paste0(mono, gsub(",", "", sub, fixed = TRUE))
  root <- which(igraph::degree(g, mode = "in") == 0L)
  fail <- function(message) {
    cli::cli_abort(message, class = "glyfp_error_unsupported")
  }
  if (
    length(root) != 1L ||
      !identical(
        igraph::graph_attr(g, "anomer"),
        paste0("?", glyrepr::infer_anomer_pos(mono[root]))
      )
  ) {
    fail("GIFFLAR requires an unspecified normal reducing end.")
  }
  if (any(!labels %in% templates[[1L]])) {
    fail(
      "GIFFLAR input contains a residue or modification outside the supported template set; see the Get started article."
    )
  }
  endpoints <- igraph::as_edgelist(g, names = FALSE)
  linkage <- igraph::edge_attr(g, "linkage")
  if (any(!grepl("^[ab][0-9]+-[0-9]+$", linkage))) {
    fail("GIFFLAR requires concrete internal linkages.")
  }
  splice <- function(x, slot, piece) {
    start <- regexpr(slot, x$text, fixed = TRUE)[1L]
    if (start < 1L) {
      fail("GIFFLAR linkage has no free attachment position.")
    }
    end <- start + nchar(slot)
    before <- if (start > 1L) seq_len(start - 1L) else integer()
    after <- if (end <= nchar(x$text)) {
      seq.int(end, nchar(x$text))
    } else {
      integer()
    }
    list(
      text = paste0(
        substr(x$text, 1L, start - 1L),
        piece$text,
        substring(x$text, end)
      ),
      owners = c(x$owners[before], piece$owners, x$owners[after])
    )
  }
  build <- function(node, anomer, depth) {
    row <- templates[
      templates[[1L]] == labels[node] & templates[[2L]] == anomer,
      ,
      drop = FALSE
    ]
    if (nrow(row) != 1L || depth > 99L) {
      fail("GIFFLAR cannot construct this residue template or ring depth.")
    }
    ring <- if (depth < 10L) as.character(depth) else sprintf("%%%02d", depth)
    text <- gsub("{ring}", ring, row[[3L]], fixed = TRUE)
    result <- list(text = text, owners = rep(node, nchar(text)))
    hetero <- strsplit(row[[4L]], ",", fixed = TRUE)[[1L]]
    hetero <- stats::setNames(sub(".*:", "", hetero), sub(":.*", "", hetero))
    allowed <- strsplit(row[[6L]], ",", fixed = TRUE)[[1L]]
    for (edge in which(endpoints[, 1L] == node)) {
      child <- endpoints[edge, 2L]
      link <- linkage[edge]
      donor <- sub("^[ab]([0-9]+)-.*", "\\1", link)
      acceptor <- sub(".*-", "", link)
      child_row <- templates[
        templates[[1L]] == labels[child] &
          templates[[2L]] == substr(link, 1L, 1L),
        ,
        drop = FALSE
      ]
      if (!identical(donor, child_row[[5L]]) || !acceptor %in% allowed) {
        fail(
          "GIFFLAR linkage uses an unsupported donor or attachment position; anomeric-to-anomeric and modification-dependent attachments are not supported."
        )
      }
      piece <- build(child, substr(link, 1L, 1L), depth + 1L)
      if (!grepl("^[ON][^=%0-9]", piece$text)) {
        fail("GIFFLAR donor has no free anomeric attachment.")
      }
      result <- splice(result, paste0("{p", acceptor, "}"), piece)
    }
    slots <- regmatches(result$text, gregexpr("\\{p[0-9]\\}", result$text))[[
      1L
    ]]
    for (slot in slots) {
      atom <- unname(hetero[substr(slot, 3L, 3L)])
      if (is.na(atom)) {
        fail("GIFFLAR template has an undefined substitution slot.")
      }
      result <- splice(result, slot, list(text = atom, owners = node))
    }
    result
  }
  result <- build(root, "?", 1L)
  result$labels <- labels
  result$linkage <- linkage
  result$endpoints <- endpoints
  result
}

gifflar_parse_smiles <- function(fragment) {
  # A restricted parser for the fixed templates, not a general user SMILES parser.
  matches <- gregexpr(
    "\\[[^]]*\\]|Br|Cl|[BCNOPSFI]|%[0-9][0-9]|[0-9]|[()=#]",
    fragment$text,
    perl = TRUE
  )[[1L]]
  tokens <- regmatches(fragment$text, list(matches))[[1L]]
  if (!identical(paste0(tokens, collapse = ""), fragment$text)) {
    cli::cli_abort(
      "GIFFLAR encountered an unsupported template token.",
      class = "glyfp_error_backend"
    )
  }
  elements <- chirality <- character()
  owners <- parent <- parent_bond <- integer()
  bonds <- matrix(integer(), ncol = 2L)
  rings <- list()
  pending <- list()
  branches <- integer()
  previous <- 0L
  for (i in seq_along(tokens)) {
    token <- tokens[i]
    if (token == "(") {
      branches <- c(branches, previous)
    } else if (token == ")") {
      previous <- utils::tail(branches, 1L)
      branches <- utils::head(branches, -1L)
    } else if (token %in% c("=", "#")) {
      # Bond orders are not features in this upstream GIFFLAR checkpoint.
      next
    } else if (grepl("^(%|[0-9])", token)) {
      if (is.null(pending[[token]])) {
        pending[[token]] <- previous
      } else {
        opened <- pending[[token]]
        pending[[token]] <- NULL
        up <- list()
        walk <- previous
        trail <- integer()
        while (walk > 0L) {
          up[[as.character(walk)]] <- trail
          trail <- c(trail, parent_bond[walk])
          walk <- parent[walk]
        }
        walk <- opened
        path <- integer()
        while (!as.character(walk) %in% names(up)) {
          path <- c(path, parent_bond[walk])
          walk <- parent[walk]
        }
        bonds <- rbind(bonds, c(opened, previous))
        rings[[length(rings) + 1L]] <- c(
          path,
          up[[as.character(walk)]],
          nrow(bonds)
        )
      }
    } else {
      element <- if (startsWith(token, "[")) {
        sub("^\\[([A-Z][a-z]?).*", "\\1", token)
      } else {
        token
      }
      stereo <- if (grepl("@@", token, fixed = TRUE)) {
        "@@"
      } else if (grepl("@", token, fixed = TRUE)) {
        "@"
      } else {
        ""
      }
      elements <- c(elements, element)
      chirality <- c(chirality, stereo)
      owners <- c(owners, fragment$owners[matches[i]])
      parent <- c(parent, previous)
      parent_bond <- c(parent_bond, if (previous > 0L) nrow(bonds) + 1L else 0L)
      if (previous > 0L) {
        bonds <- rbind(bonds, c(previous, length(elements)))
      }
      previous <- length(elements)
    }
  }
  if (length(pending) || length(branches)) {
    cli::cli_abort(
      "Invalid GIFFLAR template ring or branch.",
      class = "glyfp_error_backend"
    )
  }
  list(
    elements = elements,
    chirality = chirality,
    owners = owners,
    bonds = bonds,
    rings = rings
  )
}

gifflar_graph <- function(g, templates, vocabulary) {
  fragment <- gifflar_smiles(g, templates)
  mol <- gifflar_parse_smiles(fragment)
  ids <- function(labels) {
    value <- unname(vocabulary[labels])
    if (anyNA(value)) {
      cli::cli_abort(
        "GIFFLAR input contains a token outside the fixed vocabulary.",
        class = "glyfp_error_unknown"
      )
    }
    as.integer(value)
  }
  nb <- nrow(mol$bonds)
  atoms_to_bonds <- cbind(as.vector(t(mol$bonds)), rep(seq_len(nb), each = 2L))
  adjacent <- function(indices) {
    if (length(indices) < 2L) {
      return(matrix(integer(), ncol = 2L))
    }
    pairs <- cbind(
      rep(indices, each = length(indices)),
      rep(indices, length(indices))
    )
    pairs[pairs[, 1L] != pairs[, 2L], , drop = FALSE]
  }
  incident <- lapply(seq_along(mol$elements), function(i) {
    which(mol$bonds[, 1L] == i | mol$bonds[, 2L] == i)
  })
  atom_codes <- c(C = 1L, N = 2L, O = 3L, P = 4L, S = 5L)
  stereo_codes <- stats::setNames(1:3, c("@", "@@", ""))
  atom_x <- unname(atom_codes[mol$elements])
  bond_x <- unname(stereo_codes[match(
    mol$chirality[mol$bonds[, 1L]],
    names(stereo_codes)
  )])
  if (anyNA(c(atom_x, bond_x))) {
    cli::cli_abort(
      "GIFFLAR encountered an unsupported atom feature.",
      class = "glyfp_error_backend"
    )
  }
  # Upstream traversal writes donor nodes before acceptors. In the directed
  # glyrepr graph parents are acceptors, so reverse the mono-level edges first.
  mono_edges <- fragment$endpoints[, 2:1, drop = FALSE]
  edges <- list(
    cbind(
      as.vector(t(mol$bonds)),
      as.vector(t(mol$bonds[, 2:1, drop = FALSE]))
    ),
    atoms_to_bonds,
    cbind(seq_len(nb), mol$owners[mol$bonds[, 1L]]),
    do.call(rbind, lapply(incident, adjacent)),
    do.call(rbind, lapply(mol$rings, adjacent)),
    rbind(mono_edges, mono_edges[, 2:1, drop = FALSE])
  )
  names(edges) <- gifflar_relations()
  edges <- lapply(edges, function(x) {
    if (is.null(x)) matrix(integer(), ncol = 2L) else unname(x - 1L)
  })
  link_ids <- ids(fragment$linkage)
  list(
    x = list(atoms = atom_x, bonds = bond_x, monosacchs = ids(fragment$labels)),
    edges = edges,
    edge_attr = stats::setNames(
      list(c(link_ids, link_ids + length(vocabulary) + 1L)),
      utils::tail(gifflar_relations(), 1L)
    )
  )
}

gifflar_inputs <- function(glycans, templates, vocabulary) {
  lapply(as.list(normalize_glycans(glycans, "propagate")), function(g) {
    if (is.null(g)) {
      return(NULL)
    }
    gifflar_graph(
      as.list(glyrepr::glycan_structure(g))[[1L]],
      templates,
      vocabulary
    )
  })
}
