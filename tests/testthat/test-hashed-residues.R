# Exercise both public circular APIs in the shared mode-based checks.
circular_encode <- function(glycans, ..., mode = "binary") {
  encode <- if (mode == "binary") fp_circular_binary else fp_circular_count
  encode(glycans, ...)
}

test_that("hashed residue methods validate parameters and complete structures", {
  for (encode in list(fp_circular_binary, fp_circular_count, fp_residue_pair)) {
    for (value in list(NULL, NA, Inf, -1, 11, 1.2, "2", c(1, 2))) {
      expect_error(
        encode("Gal(b1-", radius = value),
        class = "glyfp_error_argument"
      )
    }
    for (value in list(NULL, NA, Inf, 0, 1048577, 2.5, "2048", c(1, 2))) {
      expect_error(
        encode("Gal(b1-", n_bits = value),
        class = "glyfp_error_argument"
      )
    }
    expect_error(
      encode("Gal(b1-", na_action = "omit"),
      class = "glyfp_error_argument"
    )
    expect_error(encode("Gal(b1-", motif_dict = "dynamic"), "must be empty")
    expect_error(encode("Gal(b1-", 2), "must be empty")
    for (text in c(
      "Gal",
      "Hex(b1-",
      "Gal(?1-4)GlcNAc(b1-",
      "Gal(b1-?)GlcNAc(b1-",
      "Gal?S(b1-"
    )) {
      expect_error(encode(text), class = "glyfp_error_unsupported")
      expect_error(
        encode(glyparse::auto_parse(text)),
        class = "glyfp_error_unsupported"
      )
    }
  }
})

test_that("fixed hash encodings agree with an independent SHA-256 calculation", {
  # Independently derived with Python hashlib, explicit root/child records and
  # two WL updates on Gal(b1-4)GlcNAc(b1-. No R serialization or graph API used.
  x <- "Gal(b1-4)GlcNAc(b1-"
  expected <- integer(2048)
  expected[c(345, 1836, 1828, 1191, 1280, 197)] <- 1L
  expect_identical(unname(circular_encode(x, mode = "count")[[1]]), expected)
  expected[] <- 0L
  expected[913] <- 1L
  expect_identical(unname(fp_residue_pair(x, mode = "count")[[1]]), expected)
  expect_identical(
    residue_bins(c("00000000", "ffffffff", "80000000"), 1048576L),
    c(1L, 1048576L, 1L)
  )
  expect_false(identical(
    residue_record(c("ab", "c")),
    residue_record(c("a", "bc"))
  ))
})

test_that("coordinates do not depend on batches, graph numbering or serialization", {
  x <- c(
    first = "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-",
    absent = NA_character_,
    duplicate = "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-"
  )
  for (encode in list(circular_encode, fp_residue_pair)) {
    for (mode in c("binary", "count")) {
      fp <- encode(x, mode = mode)
      expect_identical(names(fp), names(x))
      expect_identical(lengths(fp), setNames(rep(2048L, 3), names(x)))
      expect_type(fp[[1]], "integer")
      expect_identical(fp[[1]], fp[[3]])
      expect_true(all(is.na(fp[[2]])))
      expect_identical(encode(rev(x), mode = mode), rev(fp))
      expect_identical(
        c(encode(x[1], mode = mode), encode(x[2:3], mode = mode)),
        fp
      )
      objects <- glyparse::auto_parse(x)
      expect_identical(encode(objects, mode = mode), fp)
      g <- as.list(objects)[[1]]
      permuted <- glyrepr::glycan_structure(igraph::permute(
        g,
        rev(seq_len(igraph::vcount(g)))
      ))
      expect_identical(encode(permuted, mode = mode), unname(fp[1]))
      saved <- tempfile()
      on.exit(unlink(saved), add = TRUE)
      saveRDS(fp, saved)
      expect_identical(readRDS(saved), fp)
      expect_identical(encode(character(), mode = mode), list())
      expect_identical(encode(glyrepr::glycan_structure(), mode = mode), list())
      expect_identical(encode(NA_character_, mode = mode)[[1]], fp[[2]])
      expect_error(
        encode(NA_character_, na_action = "error"),
        class = "glyfp_error_missing"
      )
      # Raw helper invariance does not rely on public re-canonicalization.
      method <- if (identical(encode, circular_encode)) {
        "circular"
      } else {
        "residue_pair"
      }
      expect_identical(
        sort(residue_tokens(g, 2L, method)),
        sort(residue_tokens(
          igraph::permute(g, rev(seq_len(igraph::vcount(g)))),
          2L,
          method
        ))
      )
    }
  }
  ct <- "RES\n1b:b-dgal-HEX-1:5"
  for (encode in list(circular_encode, fp_residue_pair)) {
    expect_identical(encode(ct), encode("Gal(b1-"))
  }
})

test_that("linkages, branches, epimers, substitutions and root states are encoded", {
  x <- c(
    "Gal(b1-4)GlcNAc(b1-",
    "Gal(b1-3)GlcNAc(b1-",
    "Gal(a1-4)GlcNAc(b1-",
    "Man(b1-4)GlcNAc(b1-",
    "Gal6S(b1-4)GlcNAc(b1-",
    "Gal(b1-4)GlcNAc(a1-"
  )
  for (encode in list(circular_encode, fp_residue_pair)) {
    fp <- encode(x)
    expect_length(unique(fp), length(x))
    branches <- encode(c(
      "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-",
      "Gal(b1-6)[GlcNAc(b1-3)]GalNAc(a1-"
    ))
    expect_false(identical(branches[[1]], branches[[2]]))
  }
  # Radius-zero endpoints do not encode acceptor linkage positions.
  pair <- fp_residue_pair(x[1:2], radius = 0)
  expect_identical(pair[[1]], pair[[2]])
})

test_that("counts and binary folding have explicit occurrence and collision semantics", {
  x <- "Gal(b1-4)Gal(b1-4)Gal(b1-"
  for (radius in 0:3) {
    a <- circular_encode(x, radius = radius, mode = "count")[[1]]
    b <- fp_residue_pair(x, radius = radius, mode = "count")[[1]]
    expect_equal(sum(a), 3L * (radius + 1L))
    expect_equal(sum(b), choose(3, 2))
    expect_identical(
      unname(circular_encode(x, radius = radius)[[1]]),
      unname(as.integer(a > 0))
    )
    expect_identical(
      unname(fp_residue_pair(x, radius = radius)[[1]]),
      unname(as.integer(b > 0))
    )
    expect_equal(
      sum(circular_encode("Gal(b1-", radius = radius, mode = "count")[[1]]),
      radius + 1L
    )
  }
  for (encode in list(circular_encode, fp_residue_pair)) {
    count <- encode(x, n_bits = 1, mode = "count")[[1]]
    binary <- encode(x, n_bits = 1)[[1]]
    expect_identical(unname(binary), 1L)
    expect_equal(sum(count), if (identical(encode, circular_encode)) 9 else 3)
    # Deliberate folding collision: even different epimers become indistinguishable.
    collided <- encode(
      c("Gal(b1-4)GlcNAc(b1-", "Man(b1-4)GlcNAc(b1-"),
      n_bits = 1
    )
    expect_identical(collided[[1]], collided[[2]])
    expect_false(identical(
      names(encode(x, radius = 0)[[1]]),
      names(encode(x, radius = 1)[[1]])
    ))
  }
  expect_true(all(fp_residue_pair("Gal(b1-")[[1]] == 0L))
  expect_identical(fp_residue_pair("Gal(b1-"), fp_residue_pair("Man(a1-"))
})

test_that("controlled chains expose complementary information and limitations", {
  chain <- function(s) {
    paste0(
      paste(
        ifelse(strsplit(s, "", fixed = TRUE)[[1]] == "G", "Gal", "Man"),
        collapse = "(b1-4)"
      ),
      "(b1-"
    )
  }
  cases <- list(
    c("GGGMGGGG", "GGGGMGGG"),
    c("GGGGGGMGGGGGGGGGGGGG", "GGGGGGGMGGGGGGGGGGGG")
  )
  for (i in seq_along(cases)) {
    x <- vapply(cases[[i]], chain, "", USE.NAMES = FALSE)
    motif <- fp_motif_count(x, motif_dict = "dynamic")
    circular <- circular_encode(x, mode = "count")
    pair <- fp_residue_pair(x, mode = "count")
    expect_identical(motif[[1]], motif[[2]])
    expect_identical(identical(circular[[1]], circular[[2]]), i == 2L)
    expect_false(identical(pair[[1]], pair[[2]]))
  }
})


test_that("numeric display options do not affect feature hashes", {
  x <- "Gal(b1-4)GlcNAc(b1-4)Man(a1-"
  circular <- circular_encode(x)
  pairs <- fp_residue_pair(x)
  old <- options(scipen = -9, OutDec = ",", digits = 1)
  on.exit(options(old), add = TRUE)
  expect_identical(circular_encode(x), circular)
  expect_identical(fp_residue_pair(x), pairs)
})

test_that("mode belongs only to the residue-pair API", {
  expect_error(
    fp_residue_pair("Gal(b1-", mode = "counts"),
    class = "glyfp_error_argument"
  )
  expect_error(
    fp_residue_pair("Gal(b1-", mode = c("binary", "count")),
    class = "glyfp_error_argument"
  )
  for (encode in list(fp_circular_binary, fp_circular_count)) {
    expect_error(encode("Gal(b1-", mode = "count"), "must be empty")
  }
})
