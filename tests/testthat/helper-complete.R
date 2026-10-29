# Reference models omit the reducing-end state. Supply explicit known chemistry
# at the public boundary; only the documented model projection removes it.
complete_reference <- function(x, anomer = "b") {
  objects <- glyrepr::as_glycan_structure(x)
  graphs <- lapply(as.list(objects), function(g) {
    if (is.null(g)) {
      return(NULL)
    }
    igraph::set_graph_attr(
      g,
      "anomer",
      sub("^[?ab]", anomer, igraph::graph_attr(g, "anomer"))
    )
  })
  out <- do.call(glyrepr::glycan_structure, graphs)
  names(out) <- names(x)
  out
}

check_embedding_contract <- function(method, checkpoint, inputs, expected) {
  local_mocked_bindings(
    download.file = function(...) stop("Unexpected network"),
    .package = "utils"
  )
  before <- sweetnet_sha256(checkpoint)
  time <- file.info(checkpoint)$mtime
  x <- complete_reference(inputs)
  encode <- get(paste0("fp_", method))
  run <- function(g, path = checkpoint) {
    args <- list(glycans = g, checkpoint = path)
    if (method == "sweetnet") {
      args$download <- FALSE
    }
    do.call(encode, args)
  }
  fp <- run(x)
  expect_lt(max(abs(do.call(rbind, fp) - expected)), 1e-5)
  expect_identical(run(as.character(x)), fp)
  expect_identical(run(complete_reference(inputs, "a")), fp)
  selected <- x[c(2, 1, 2, NA_integer_)]
  names(selected) <- c("second", "first", "duplicate", "missing")
  output <- run(selected)
  expect_identical(names(output), names(selected))
  expect_identical(output[[1]], output[[3]])
  expect_identical(output[[2]], fp[[1]])
  expect_true(all(is.na(output[[4]])))
  expect_identical(run(rev(selected)), rev(output))
  expect_identical(c(run(selected[1]), run(selected[2:4])), output)
  expect_identical(run(character()), list())
  expect_identical(run(NA_character_)[[1]], output[[4]])
  saved <- tempfile()
  corrupt <- tempfile(fileext = ".pt")
  on.exit(unlink(c(saved, corrupt)), add = TRUE)
  saveRDS(output, saved)
  expect_identical(readRDS(saved), output)
  writeLines("not weights", corrupt)
  expect_error(
    run(x[1], path = corrupt),
    class = if (method == "sweetnet") {
      "glyfp_error_checksum"
    } else {
      "glyfp_error_backend"
    }
  )
  expect_identical(sweetnet_sha256(checkpoint), before)
  expect_identical(file.info(checkpoint)$mtime, time)
  set.seed(725)
  rng <- get(".Random.seed", envir = globalenv())
  torch_rng <- torch::torch_get_rng_state()$clone()
  threads <- torch::torch_get_num_threads()
  expect_identical(run(x[1]), fp[1])
  expect_identical(get(".Random.seed", envir = globalenv()), rng)
  expect_true(torch::torch_equal(torch::torch_get_rng_state(), torch_rng))
  expect_identical(torch::torch_get_num_threads(), threads)
}
