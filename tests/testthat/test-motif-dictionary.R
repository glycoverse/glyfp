test_that("motif_dict is validated without changing fixed defaults", {
  x <- c(a = "Gal(b1-3)Gal(b1-", absent = NA_character_)
  for (encode in list(fp_motif_count, fp_motif_binary)) {
    expect_identical(formals(encode)$motif_dict, "fixed")
    expect_identical(encode(x), encode(x, motif_dict = "fixed"))
    expect_length(encode(x)[[1]], 3173L)
    for (value in list(
      NULL,
      character(),
      NA_character_,
      "",
      "auto",
      "FIXED",
      "dyn",
      c("fixed", "dynamic"),
      1,
      TRUE
    )) {
      expect_error(
        encode(x, motif_dict = value),
        class = "glyfp_error_argument"
      )
    }
  }
  for (encode in list(fp_sweetnet, fp_glycangt, fp_gifflar)) {
    expect_false("motif_dict" %in% names(formals(encode)))
  }
})

test_that("dynamic counts include novel motifs and preserve overlap semantics", {
  x <- c(
    chain = "Gal(b1-3)Gal(b1-3)Gal(b1-",
    novel = "Gal(b1-9)GlcNAc(b1-",
    absent = NA_character_
  )
  fixed <- fp_motif_count(x)
  fp <- fp_motif_count(x, motif_dict = "dynamic")
  features <- names(fp[[1]])
  expect_true(x[["novel"]] %in% features)
  expect_false(x[["novel"]] %in% names(fixed[[1]]))
  expect_identical(fp[[2]][[x[["novel"]]]], 1L)
  expect_identical(fp[[1]][["Gal(b1-3)Gal(b1-"]], 2L)
  expect_identical(fp[[1]][[x[["chain"]]]], 1L)
  expected <- glymotif::count_motifs(
    glyparse::auto_parse(x[1:2]),
    glyrepr::as_glycan_structure(features),
    alignments = "substructure",
    ignore_linkages = FALSE,
    strict_sub = TRUE,
    mode = "strict"
  )
  expect_equal(unname(do.call(rbind, fp[1:2])), unname(expected))
  expect_identical(names(fp), names(x))
  expect_type(fp[[1]], "integer")
  expect_true(all(is.na(fp[[3]])))
  expect_identical(
    fp_motif_binary(x, motif_dict = "dynamic"),
    lapply(fp, \(x) 1L * (x > 0L))
  )
  for (unknown in c("ignore", "warn", "error")) {
    expect_identical(
      fp_motif_count(x, motif_dict = "dynamic", unknown = unknown),
      fp
    )
  }
  # Dynamic calls neither read nor mutate the bundled dictionary.
  expect_identical(fp_motif_count(x), fixed)
  local_mocked_bindings(motif_features = function(...) {
    stop("fixed dictionary accessed")
  })
  expect_identical(fp_motif_count(x, motif_dict = "dynamic"), fp)
})

test_that("dynamic dictionaries are canonical and deterministic within motif sets", {
  x <- c("Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-", "Man(a1-")
  for (encode in list(fp_motif_count, fp_motif_binary)) {
    fp <- encode(x, motif_dict = "dynamic")
    expect_identical(
      names(fp[[1]]),
      sort(unique(enc2utf8(names(fp[[1]]))), method = "radix")
    )
    expect_identical(
      encode(x[c(2, 1, 2)], motif_dict = "dynamic"),
      fp[c(2, 1, 2)]
    )
    objects <- glyparse::auto_parse(x)
    expect_identical(encode(objects, motif_dict = "dynamic"), fp)
    g <- as.list(objects)[[1]]
    permuted <- glyrepr::glycan_structure(igraph::permute(
      g,
      rev(seq_len(igraph::vcount(g)))
    ))
    expect_identical(
      encode(objects[1], motif_dict = "dynamic"),
      encode(permuted, motif_dict = "dynamic")
    )
    a <- encode("Gal(b1-", motif_dict = "dynamic")[[1]]
    b <- encode("Man(a1-", motif_dict = "dynamic")[[1]]
    expect_length(a, length(b))
    expect_false(identical(names(a), names(b)))
    expect_identical(
      names(encode("Gal(b1-3)Gal(b1-", motif_dict = "dynamic", max_size = 1)[[
        1
      ]]),
      names(a)
    )
    expect_identical(
      names(encode("Gal(b1-3)Gal(b1-", motif_dict = "dynamic", max_size = 2)[[
        1
      ]]),
      c("Gal(b1-", "Gal(b1-3)Gal(b1-")
    )
    ends <- encode(c("Gal(a1-", "Gal(b1-"), motif_dict = "dynamic")
    expect_false(identical(ends[[1]], ends[[2]]))
    expect_identical(ends[[1]], c("Gal(a1-" = 1L, "Gal(b1-" = 0L))
  }
})

test_that("dynamic empty and missing inputs have explicit zero-dimensional behavior", {
  for (encode in list(fp_motif_count, fp_motif_binary)) {
    for (x in list(character(), glyrepr::glycan_structure())) {
      expect_identical(encode(x, motif_dict = "dynamic"), list())
    }
    missing <- c(first = NA_character_, second = NA_character_)
    fp <- encode(missing, motif_dict = "dynamic")
    expect_identical(names(fp), names(missing))
    expect_identical(
      fp,
      setNames(rep(list(setNames(integer(), character())), 2), names(missing))
    )
    mixed <- encode(c(missing, "Gal(b1-"), motif_dict = "dynamic")
    expect_identical(mixed[[1]], c("Gal(b1-" = NA_integer_))
    expect_error(
      encode(missing, motif_dict = "dynamic", na_action = "error"),
      class = "glyfp_error_missing"
    )
    for (x in c("Gal", "Hex(b1-", "Gal(b1-?)GlcNAc(b1-")) {
      expect_error(
        encode(x, motif_dict = "dynamic"),
        class = "glyfp_error_unsupported"
      )
      expect_error(
        encode(glyparse::auto_parse(x), motif_dict = "dynamic"),
        class = "glyfp_error_unsupported"
      )
    }
  }
})
