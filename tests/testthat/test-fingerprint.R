test_that("the sole public API uses fixed coordinates across independent batches", {
  expect_identical(getNamespaceExports("glyfp"), "glycan_fingerprint")
  x <- c(
    first = "Gal(b1-3)Gal(b1-3)Gal(b1-",
    absent = NA,
    duplicate = "Galb3Galb3Galb-"
  )
  fp <- glycan_fingerprint(x, "substructure_count")
  expect_type(fp, "list")
  expect_identical(names(fp), names(x))
  expect_identical(lengths(fp), setNames(rep(3173L, 3), names(x)))
  expect_identical(fp[[1]], fp[[3]])
  expect_true(all(is.na(fp[[2]])))
  expect_identical(
    fp[[1]][c("Gal(b1-", "Gal(b1-3)Gal(b1-")],
    c("Gal(b1-" = 3L, "Gal(b1-3)Gal(b1-" = 2L)
  )
  expect_identical(
    c(
      glycan_fingerprint(x[1], "substructure_count"),
      glycan_fingerprint(x[2:3], "substructure_count")
    ),
    fp
  )
  expect_identical(glycan_fingerprint(rev(x), "substructure_count"), rev(fp))
  expect_identical(names(glycan_fingerprint("Man(a1-")[[1]]), names(fp[[1]]))
  expect_identical(glycan_fingerprint(x), lapply(fp, function(x) 1L * (x > 0L)))
})

test_that("candidate pruning preserves exhaustive glymotif matching and overlaps", {
  x <- c(
    "Gal(b1-3)Gal(b1-3)Gal(b1-",
    "Gal(b1-3)[Gal(b1-6)]GlcNAc(b1-",
    "Gal6S(b1-4)GlcNAc(a1-"
  )
  actual <- glycan_fingerprint(x, "substructure_count")
  expected <- glymotif::count_motifs(
    glyparse::auto_parse(x),
    glyrepr::as_glycan_structure(names(actual[[1]])),
    alignments = "substructure",
    ignore_linkages = FALSE,
    strict_sub = TRUE,
    mode = "strict"
  )
  expect_equal(unname(do.call(rbind, actual)), unname(expected))
  expect_identical(actual[[1]][["Gal(b1-3)Gal(b1-"]], 2L)
  # Identical branches exercise automorphism deduplication at the private engine.
  symmetric <- glyparse::auto_parse("Gal(b1-?)[Gal(b1-?)]GlcNAc(b1-")
  expect_identical(
    as.integer(glymotif::count_motifs(
      symmetric,
      symmetric,
      alignments = "substructure",
      ignore_linkages = FALSE,
      strict_sub = TRUE,
      mode = "strict"
    )),
    1L
  )
})

test_that("dictionary provenance and full feature order are stable", {
  resource <- system.file("extdata", "motif-dictionary.tsv", package = "glyfp")
  metadata <- dget(system.file(
    "extdata",
    "motif-dictionary-metadata.R",
    package = "glyfp"
  ))
  expect_identical(metadata$id, "glydb-2.12.1-complete-v1")
  expect_identical(
    metadata$included_accessions + metadata$excluded_accessions,
    metadata$source_rows
  )
  expect_identical(metadata$unique_structures, 8211L)
  expect_identical(
    digest::digest(file = resource, algo = "sha256"),
    metadata$dictionary_sha256
  )
  features <- motif_features(7L)
  expect_length(features, 24124L)
  expect_identical(features, sort(unique(enc2utf8(features)), method = "radix"))
  for (size in 1:7) {
    expect_length(
      motif_features(size),
      sum(unlist(metadata$sizes)[seq_len(size)])
    )
  }
  expect_length(glycan_fingerprint("Gal(b1-", max_size = 1)[[1]], 134L)
  expect_identical(glycan_fingerprint("Gal(b1-")[[1]][["Gal(b1-"]], 1L)
})

test_that("complete structure validation precedes every model and download", {
  local_mocked_bindings(
    sweetnet_run = function(...) stop("backend must not run"),
    glycangt_backend = function(...) stop("backend must not run"),
    gifflar_backend = function(...) stop("backend must not run")
  )
  uncertain <- c(
    "Gal",
    "Gal(?1-",
    "Gal(??-",
    "Gal(b1-?)GlcNAc(b1-",
    "Gal(?1-4)GlcNAc(b1-",
    "Hex(b1-",
    "Gal?S(b1-",
    "{Neu5Ac(a2-3)|2,3}Gal(b1-3)[Gal(b1-4)]GlcNAc(a1-"
  )
  for (text in uncertain) {
    for (method in c(
      "substructure_binary",
      "substructure_count",
      "sweetnet",
      "glycangt",
      "gifflar"
    )) {
      args <- list(method = method)
      if (method %in% c("glycangt", "gifflar")) {
        args$checkpoint <- "unused.pt"
      }
      for (x in list(text, glyparse::auto_parse(text))) {
        expect_error(
          do.call(glycan_fingerprint, c(list(glycans = x), args)),
          class = "glyfp_error_unsupported",
          info = paste(text, method)
        )
      }
    }
  }
})

test_that("reducing end chemistry is retained in motif fingerprints", {
  fp <- glycan_fingerprint(c("Gal(a1-", "Gal(b1-"))
  expect_false(identical(fp[[1]], fp[[2]]))
  expect_identical(
    fp[[1]][c("Gal(a1-", "Gal(b1-")],
    c("Gal(a1-" = 1L, "Gal(b1-" = 0L)
  )
  expect_identical(
    fp[[2]][c("Gal(a1-", "Gal(b1-")],
    c("Gal(a1-" = 0L, "Gal(b1-" = 1L)
  )
  x <- glyparse::auto_parse(c("Gal(b1-4)GlcNAc(a1-", "Gal(b1-4)GlcNAc(b1-"))
  expect_identical(
    as.character(embedding_projection(x))[1],
    as.character(embedding_projection(x))[2]
  )
  expect_identical(
    as.character(x),
    c("Gal(b1-4)GlcNAc(a1-", "Gal(b1-4)GlcNAc(b1-")
  )
})

test_that("unseen complete motifs never expand the dictionary", {
  x <- "Gal(b1-9)GlcNAc(b1-"
  expect_false(x %in% motif_features(3L))
  fp <- glycan_fingerprint(x)
  expect_length(fp[[1]], 3173L)
  expect_identical(fp[[1]][["Gal(b1-"]], 1L)
  expect_warning(
    warned <- glycan_fingerprint(x, unknown = "warn"),
    class = "glyfp_warning_unknown"
  )
  expect_identical(warned, fp)
  expect_error(
    glycan_fingerprint(x, unknown = "error"),
    class = "glyfp_error_unknown"
  )
})

test_that("text formats, graph numbering, missingness and serialization are stable", {
  x <- c(one = "Galb3GalNAca-", two = "Gal(b1-4)GlcNAc(b1-")
  object <- glyparse::auto_parse(x)
  expect_identical(glycan_fingerprint(x), glycan_fingerprint(object))
  g <- as.list(object)[[1]]
  permuted <- glyrepr::glycan_structure(igraph::permute(
    g,
    rev(seq_len(igraph::vcount(g)))
  ))
  expect_identical(
    unname(glycan_fingerprint(object[1])),
    glycan_fingerprint(permuted)
  )
  for (text in c(
    "RES\n1b:b-dgal-HEX-1:5",
    "WURCS=2.0/2,3,2/[a2112h-1b_1-5][a2112h-1a_1-4]/1-2-2/a4-b1_b2-c1"
  )) {
    expect_identical(
      glycan_fingerprint(text),
      glycan_fingerprint(glyparse::auto_parse(text))
    )
  }
  expect_identical(glycan_fingerprint(character()), list())
  expect_identical(glycan_fingerprint(glyrepr::glycan_structure()), list())
  expect_identical(
    glycan_fingerprint(NA_character_)[[1]],
    setNames(rep(NA_integer_, 3173), motif_features(3L))
  )
  expect_error(
    glycan_fingerprint(NA_character_, na_action = "error"),
    class = "glyfp_error_missing"
  )
  set.seed(17)
  rng <- .Random.seed
  fp <- glycan_fingerprint(x)
  expect_identical(glycan_fingerprint(x), fp)
  expect_identical(.Random.seed, rng)
  path <- tempfile()
  on.exit(unlink(path))
  saveRDS(fp, path)
  expect_identical(readRDS(path), fp)
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
    expect_error(glycan_fingerprint(x), class = "glyfp_error_unsupported")
  }
  bad <- glyrepr::as_glycan_structure("Gal(b1-")
  attr(bad, "graphs") <- list()
  expect_error(glycan_fingerprint(bad), class = "glyfp_error_input")
})

test_that("invalid method options cannot be silently ignored", {
  for (value in list(0, -1, 1.2, 8, Inf, NA, "3", c(1, 2))) {
    expect_error(
      glycan_fingerprint("Gal", max_size = value),
      class = "glyfp_error_argument"
    )
  }
  expect_error(
    glycan_fingerprint("Gal", method = "embedding"),
    class = "glyfp_error_argument"
  )
  expect_error(
    glycan_fingerprint("Gal", unknown = "grow"),
    class = "glyfp_error_argument"
  )
  expect_error(
    glycan_fingerprint("Gal", na_action = "drop"),
    class = "glyfp_error_argument"
  )
  expect_error(glycan_fingerprint("Gal", max_siz = 2), "must be empty")
  expect_error(
    glycan_fingerprint("Gal", "substructure_binary", 2),
    "must be empty"
  )
})
