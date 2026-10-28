test_that("API returns only a list preserving order, names and missingness", {
  x <- c(
    first = "Gal(b1-3)Gal(b1-3)Gal(b1-",
    absent = NA,
    duplicate = "Galb3Galb3Galb-"
  )
  e <- fit_fingerprint(x, "substructure_count", max_size = 2)
  fp <- glycan_fingerprint(x, "substructure_count", max_size = 2)
  expect_type(fp, "list")
  expect_identical(names(fp), names(x))
  expect_length(fp, 3)
  expect_identical(fp[[1]], c("Gal(b1-" = 3L, "Gal(b1-3)Gal(b1-" = 2L))
  expect_identical(fp[[1]], fp[[3]])
  expect_true(all(is.na(fp[[2]])))
  expect_identical(fp, transform_fingerprint(x, e))
  expect_identical(
    as.list(formals(transform_fingerprint)),
    alist(glycans = , encoder = )
  )
  expect_true(all(vapply(
    fp,
    function(x) identical(names(x), e$features),
    TRUE
  )))
})

test_that("the default includes every size up to three, including short glycans", {
  x <- c("Gal(b1-", "Gal(b1-3)Gal(b1-", "Gal(b1-3)Gal(b1-3)Gal(b1-")
  fp <- glycan_fingerprint(x, "substructure_count")
  features <- c("Gal(b1-", "Gal(b1-3)Gal(b1-", "Gal(b1-3)Gal(b1-3)Gal(b1-")
  expect_identical(fp[[1]], stats::setNames(c(1L, 0L, 0L), features))
  expect_identical(fp[[2]], stats::setNames(c(2L, 1L, 0L), features))
  expect_identical(fp[[3]], stats::setNames(c(3L, 2L, 1L), features))
})

test_that("overlap, branches and symmetry use glymotif counts", {
  x <- c("Gal(b1-3)Gal(b1-3)Gal(b1-", "Gal(b1-?)[Gal(b1-?)]GlcNAc(b1-")
  e <- fit_fingerprint(x, "substructure_count")
  counts <- transform_fingerprint(x, e)
  expected <- glymotif::count_motifs(
    glyparse::auto_parse(x),
    glyrepr::as_glycan_structure(e$features),
    alignments = "substructure",
    ignore_linkages = FALSE,
    strict_sub = TRUE,
    mode = "strict"
  )
  expect_equal(unname(do.call(rbind, counts)), unname(expected))
  expect_identical(
    glycan_fingerprint(x),
    lapply(counts, function(x) 1L * (x > 0L))
  )
  expect_identical(counts[[1]][["Gal(b1-3)Gal(b1-"]], 2L)
  expect_identical(counts[[2]][["Gal(b1-?)[Gal(b1-?)]GlcNAc(b1-"]], 1L)
  expect_identical(counts[[2]][["Gal(b1-?)GlcNAc(b1-"]], 2L)
})

test_that("dictionaries are invariant to input order and graph numbering", {
  x <- glyparse::auto_parse(c("Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-", "Man(a1-"))
  e <- fit_fingerprint(x)
  expect_identical(e, fit_fingerprint(x[c(2, 1, 2)]))
  g <- as.list(x)[[1]]
  permuted <- glyrepr::glycan_structure(igraph::permute(
    g,
    rev(seq_len(igraph::vcount(g)))
  ))
  expect_identical(
    transform_fingerprint(x[1], e),
    transform_fingerprint(permuted, e)
  )
  expect_identical(e$features, sort(e$features, method = "radix"))
})

test_that("text formats and objects share an encoder", {
  text <- c(one = "Galb3GalNAca-", two = "Gal(b1-4)GlcNAc")
  structures <- glyparse::auto_parse(text)
  e <- fit_fingerprint(structures)
  expect_identical(
    transform_fingerprint(text, e),
    transform_fingerprint(structures, e)
  )
  expect_identical(fit_fingerprint(text), e)
  expect_identical(glycan_fingerprint(text), glycan_fingerprint(structures))
  ct <- "RES\n1b:b-dgal-HEX-1:5"
  expect_identical(
    glycan_fingerprint(ct),
    glycan_fingerprint(glyparse::parse_glycoct(ct))
  )
  expect_error(glycan_fingerprint("not-a-glycan"), class = "glyfp_error_parse")
  expect_error(
    transform_fingerprint(c("Galb-", "bad input"), e),
    class = "glyfp_error_parse"
  )
  for (x in list(1:2, matrix("Gal", 1), NULL)) {
    expect_error(fit_fingerprint(x), class = "glyfp_error_input")
  }
})

test_that("equal lengths from independent fits do not imply common features", {
  a <- "Gal(b1-"
  b <- "Man(a1-"
  expect_identical(
    unname(glycan_fingerprint(a)[[1]]),
    unname(glycan_fingerprint(b)[[1]])
  )
  expect_false(identical(
    names(glycan_fingerprint(a)[[1]]),
    names(glycan_fingerprint(b)[[1]])
  ))
  e <- fit_fingerprint(c(a, b))
  expect_identical(
    c(transform_fingerprint(a, e), transform_fingerprint(b, e)),
    glycan_fingerprint(c(a, b))
  )
  expect_identical(
    transform_fingerprint(b, fit_fingerprint(a))[[1]],
    c("Gal(b1-" = 0L)
  )
})

test_that("unknown policies leave the dictionary unchanged", {
  for (policy in c("ignore", "warn", "error")) {
    e <- fit_fingerprint("Gal(b1-", unknown = policy)
    before <- serialize(e, NULL)
    if (policy == "ignore") {
      expect_identical(
        transform_fingerprint("Man(a1-", e)[[1]],
        c("Gal(b1-" = 0L)
      )
    } else if (policy == "warn") {
      expect_warning(
        transform_fingerprint(c("Man(a1-", "Glc(a1-"), e),
        class = "glyfp_warning_unknown"
      )
    } else {
      expect_error(
        transform_fingerprint("Man(a1-", e),
        class = "glyfp_error_unknown"
      )
    }
    expect_identical(serialize(e, NULL), before)
  }
  e <- fit_fingerprint("Hex(?1-", unknown = "warn")
  expect_warning(
    out <- transform_fingerprint("Gal(b1-", e),
    class = "glyfp_warning_unknown"
  )
  expect_identical(unname(out[[1]]), 1L)
})

test_that("empty and missing inputs are explicit", {
  for (x in list(character(), glyrepr::glycan_structure())) {
    e <- fit_fingerprint(x)
    expect_identical(e$features, character())
    expect_identical(glycan_fingerprint(x), list())
    expect_identical(transform_fingerprint(x, e), list())
    expect_identical(
      transform_fingerprint("Gal(b1-", e),
      list(stats::setNames(integer(), character()))
    )
  }
  e <- fit_fingerprint(c(NA_character_, NA_character_))
  expect_identical(e$features, character())
  expect_length(glycan_fingerprint(c(NA_character_, NA_character_)), 2)
  expect_identical(
    lengths(transform_fingerprint(c(NA_character_, "Gal(b1-"), e)),
    c(0L, 0L)
  )
  e <- fit_fingerprint("Gal(b1-", na_action = "error")
  expect_error(
    transform_fingerprint(NA_character_, e),
    class = "glyfp_error_missing"
  )
  expect_error(
    fit_fingerprint(c("Gal(b1-", NA), na_action = "error"),
    class = "glyfp_error_missing"
  )
})

test_that("linkages, half-linkages and substituents follow glymotif semantics", {
  x <- c(
    "Gal(b1-3)GalNAc(a1-",
    "Gal(b1-4)GalNAc(a1-",
    "Gal(a1-3)GalNAc(a1-",
    "Gal6S(b1-3)GalNAc(a1-"
  )
  fp <- glycan_fingerprint(x)
  for (i in 2:4) {
    expect_false(identical(fp[[1]], fp[[i]]))
  }
  e <- fit_fingerprint("Gal(b1-?)GalNAc(?1-")
  expect_identical(
    transform_fingerprint(x[1], e)[[1]][["Gal(b1-?)GalNAc(?1-"]],
    1L
  )
  specific <- fit_fingerprint(x[1])
  expect_identical(
    transform_fingerprint("Gal(b1-?)GalNAc(a1-", specific)[[1]][[x[1]]],
    0L
  )
})

test_that("unsupported graphs and corrupted storage fail clearly", {
  floating <- glyrepr::as_glycan_structure(
    "{Neu5Ac(a2-3)|2,3}Gal(b1-3)[Gal(b1-4)]GlcNAc(a1-"
  )
  floating_sub <- glyrepr::as_glycan_structure(
    "{6S|1,2}Gal(a1-3)Glc(a1-3)Man(a1-"
  )
  g <- as.list(glyrepr::as_glycan_structure("Gal(b1-"))[[1]]
  g <- igraph::set_graph_attr(g, "alditol", TRUE)
  alditol <- glyrepr::glycan_structure(g)
  for (x in list(floating, floating_sub, alditol)) {
    expect_error(fit_fingerprint(x), class = "glyfp_error_unsupported")
    expect_error(
      transform_fingerprint(x, fit_fingerprint("Gal(b1-")),
      class = "glyfp_error_unsupported"
    )
  }
  bad <- glyrepr::as_glycan_structure("Gal(b1-")
  attr(bad, "graphs") <- list()
  expect_error(fit_fingerprint(bad), class = "glyfp_error_input")
})

test_that("serialization, determinism and encoder validation hold", {
  x <- c("Gal(b1-3)Gal(b1-", "Man(a1-")
  e <- fit_fingerprint(x, "substructure_count")
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(e, path)
  expect_identical(readRDS(path), e)
  expect_identical(
    transform_fingerprint(x, readRDS(path)),
    transform_fingerprint(x, e)
  )
  set.seed(17)
  seed <- .Random.seed
  expect_identical(glycan_fingerprint(x), glycan_fingerprint(x))
  expect_identical(.Random.seed, seed)
  expect_output(print(e), "glyfp_encoder")
  for (field in c(
    "features",
    "method",
    "params",
    "schema_version",
    "versions"
  )) {
    broken <- e
    broken[field] <- list(NULL)
    expect_error(
      transform_fingerprint(x, broken),
      class = "glyfp_error_encoder"
    )
  }
  broken <- e
  broken$features <- rep(broken$features[1], 2)
  expect_error(transform_fingerprint(x, broken), class = "glyfp_error_encoder")
  broken <- e
  broken$versions[1] <- "0.0.0"
  expect_error(transform_fingerprint(x, broken), class = "glyfp_error_version")
  expect_error(transform_fingerprint(x, list()), class = "glyfp_error_encoder")
})

test_that("invalid method options cannot be silently ignored", {
  for (value in list(0, -1, 1.2, 8, Inf, NA, "3", c(1, 2))) {
    expect_error(
      fit_fingerprint("Gal", max_size = value),
      class = "glyfp_error_argument"
    )
  }
  expect_error(
    fit_fingerprint("Gal", method = "embedding"),
    class = "glyfp_error_argument"
  )
  expect_error(
    fit_fingerprint("Gal", unknown = "grow"),
    class = "glyfp_error_argument"
  )
  expect_error(
    fit_fingerprint("Gal", na_action = "drop"),
    class = "glyfp_error_argument"
  )
  expect_error(fit_fingerprint("Gal", max_siz = 2), "must be empty")
  expect_error(
    glycan_fingerprint("Gal", "substructure_binary", 2),
    "must be empty"
  )
})

test_that("WURCS input shares canonical coordinates", {
  wurcs <- paste0(
    "WURCS=2.0/2,3,2/",
    "[a2112h-1b_1-5][a2112h-1a_1-4]/",
    "1-2-2/a4-b1_b2-c1"
  )
  object <- glyparse::parse_wurcs(wurcs)
  e <- fit_fingerprint(object, "substructure_count")
  expect_identical(
    transform_fingerprint(wurcs, e),
    transform_fingerprint(object, e)
  )
  expect_identical(
    glycan_fingerprint(wurcs),
    glycan_fingerprint(as.character(object))
  )
})

test_that("max_size controls the feature dictionary", {
  x <- "Gal(b1-3)Gal(b1-3)Gal(b1-"
  expect_identical(fit_fingerprint(x, max_size = 1)$features, "Gal(b1-")
  expect_length(fit_fingerprint(x, max_size = 2)$features, 2)
  expect_identical(
    fit_fingerprint(x, max_size = 7)$features,
    fit_fingerprint(x, max_size = 3)$features
  )
})
