test_that("motif functions validate their own options", {
  for (encode in list(fp_motif_binary, fp_motif_count)) {
    for (value in list(0, -1, 1.2, 8, Inf, NA, "3", c(1, 2))) {
      expect_error(
        encode("Gal", max_size = value),
        class = "glyfp_error_argument"
      )
    }
    expect_error(
      encode("Gal", unknown = "grow"),
      class = "glyfp_error_argument"
    )
    expect_error(
      encode("Gal", na_action = "drop"),
      class = "glyfp_error_argument"
    )
    expect_error(
      encode(NA_character_, na_action = "error"),
      class = "glyfp_error_missing"
    )
    expect_identical(encode(character()), list())
    expect_length(encode(NA_character_, max_size = 1)[[1]], 134L)
    expect_error(encode("Gal", max_siz = 2), "must be empty")
    expect_error(encode("Gal", method = "substructure_count"), "must be empty")
    expect_error(encode("Gal", 2), "must be empty")
    expect_error(
      encode("Gal(b1-9)GlcNAc(b1-", unknown = "error"),
      class = "glyfp_error_unknown"
    )
  }
})

test_that("embedding APIs preserve projection, batches, names and missingness", {
  # Exercise each public adapter without downloading weights or loading torch.
  seen <- list()
  produce <- function(x, n) {
    seen[[length(seen) + 1L]] <<- as.character(x)
    lapply(as.character(x), \(key) seq_len(n) / 1000 + sum(utf8ToInt(key)))
  }
  local_mocked_bindings(
    sweetnet_run = function(inputs, checkpoint, cache_dir, download) {
      expect_identical(checkpoint, "local.pt")
      produce(vapply(inputs, `[[`, "", "iupac"), 128L)
    },
    glycangt_backend = function(glycans, checkpoint) {
      expect_identical(checkpoint, "local.pt")
      produce(glycans, 768L)
    },
    gifflar_backend = function(glycans, checkpoint) {
      expect_identical(checkpoint, "local.pt")
      produce(glycans, 384L)
    }
  )
  x <- c(
    first = "Gal(b1-4)GlcNAc(b1-",
    second = "Gal(b1-3)GlcNAc(b1-",
    same_root_projection = "Gal(b1-4)GlcNAc(a1-",
    missing = NA_character_,
    duplicate = "Gal(b1-4)GlcNAc(b1-"
  )
  for (method in c("sweetnet", "glycangt", "gifflar")) {
    encode <- get(paste0("fp_", method))
    run <- function(x) encode(x, checkpoint = "local.pt")
    fp <- run(x)
    projected <- seen[[length(seen)]]
    expect_length(projected, 2L)
    expected <- embedding_projection(glyparse::auto_parse(x[1:2]))
    if (method == "sweetnet") {
      expected <- vapply(sweetnet_inputs(expected), `[[`, "", "iupac")
    }
    expect_identical(unname(projected), unname(as.character(expected)))
    expect_identical(names(fp), names(x))
    features <- get(paste0(method, "_features"))()
    expect_identical(names(fp[[1]]), features)
    expect_type(fp[[1]], "double")
    expect_identical(fp[[1]], fp[[3]])
    expect_identical(fp[[1]], fp[[5]])
    expect_false(identical(fp[[1]], fp[[2]]))
    expect_identical(
      fp[[4]],
      setNames(rep(NA_real_, length(features)), features)
    )
    expect_identical(run(rev(x)), rev(fp))
    expect_identical(c(run(x[1]), run(x[-1])), fp)
    expect_identical(run(glyparse::auto_parse(x)), fp)
    expect_identical(run(character()), list())
    expect_identical(run(NA_character_)[[1]], fp[[4]])
  }
})

test_that("SweetNet owns checkpoint, cache and download defaults", {
  calls <- list()
  local_mocked_bindings(sweetnet_run = function(
    inputs,
    checkpoint,
    cache_dir,
    download
  ) {
    calls[[length(calls) + 1L]] <<- list(checkpoint, cache_dir, download)
    lapply(inputs, \(x) rep(0, 128L))
  })
  fp_sweetnet("Gal(b1-")
  expect_identical(
    calls[[1]],
    list("sweetnet-v2", tools::R_user_dir("glyfp", "cache"), TRUE)
  )
  fp_sweetnet(
    "Gal(b1-",
    checkpoint = "sweetnet-v1",
    cache_dir = "custom",
    download = FALSE
  )
  expect_identical(calls[[2]], list("sweetnet-v1", "custom", FALSE))
})

test_that("embedding adapters reject invalid backend results", {
  output <- NULL
  local_mocked_bindings(sweetnet_run = function(...) output)
  for (bad in list(
    list(),
    list(rep(0, 127)),
    list(rep(NA_real_, 128)),
    list(rep(Inf, 128)),
    list(rep(0, 128), rep(0, 128))
  )) {
    output <- bad
    expect_error(fp_sweetnet("Gal(b1-"), class = "glyfp_error_backend")
  }
})
