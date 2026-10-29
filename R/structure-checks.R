# Inspect normalized object fields, not the input string. Empty substituent and
# floating-component fields mean absence; they are not unknown chemistry.
structure_issues <- function(graph) {
  reasons <- character()
  if (glyrepr::get_mono_type(graph) != "concrete") {
    reasons <- c(reasons, "nonconcrete_residue")
  }
  if (length(igraph::graph_attr(graph, "floating_parts"))) {
    reasons <- c(reasons, "floating_parts")
  }
  if (length(igraph::graph_attr(graph, "floating_substituents"))) {
    reasons <- c(reasons, "floating_substituents")
  }
  alditol <- isTRUE(igraph::graph_attr(graph, "alditol"))
  if (alditol) {
    reasons <- c(reasons, "unsupported_alditol")
  }
  anomer <- igraph::graph_attr(graph, "anomer")
  # Alditols have no anomeric stereocentre, but remain outside package scope.
  if (
    !alditol &&
      (length(anomer) != 1L ||
        is.na(anomer) ||
        !grepl("^[ab][1-9][0-9]*$", anomer))
  ) {
    reasons <- c(reasons, "uncertain_reducing_end")
  }
  links <- igraph::edge_attr(graph, "linkage")
  if (anyNA(links) || any(!grepl("^[ab][1-9][0-9]*-[1-9][0-9]*$", links))) {
    reasons <- c(reasons, "uncertain_linkage")
  }
  subs <- igraph::vertex_attr(graph, "sub")
  if (anyNA(subs) || any(grepl("[?/]", subs))) {
    reasons <- c(reasons, "uncertain_substituent")
  }
  reasons
}

# User-approved model projection: the input has already passed complete-structure
# validation. Only the known reducing-end alpha/beta state is omitted, because
# the pretrained models do not encode it. Preserve its carbon position, allowing
# each backend to reject non-normal reducing-end positions as before.
embedding_projection <- function(glycans) {
  if (!length(glycans)) {
    return(glycans)
  }
  graphs <- lapply(as.list(glycans), function(graph) {
    if (is.null(graph)) {
      return(NULL)
    }
    igraph::set_graph_attr(
      graph,
      "anomer",
      sub("^[ab]", "?", igraph::graph_attr(graph, "anomer"))
    )
  })
  result <- do.call(glyrepr::glycan_structure, graphs)
  names(result) <- names(glycans)
  result
}
