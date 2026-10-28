test_that("SweetNet normalization only removes unknown inferred half-linkages", {
  x <- glyparse::auto_parse(c("Gal(b1-4)GlcNAc", "GlcNAc", "Neu5Ac"))
  inputs <- sweetnet_inputs(x)
  expect_identical(
    vapply(inputs, `[[`, "", "iupac"),
    c("Gal(b1-4)GlcNAc", "GlcNAc", "Neu5Ac")
  )
  expect_identical(
    sweetnet_inputs(glyparse::auto_parse(as.character(x))),
    inputs
  )
  expect_identical(inputs[[2]]$edges, list())
  expect_identical(as.character(inputs[[1]]$labels), c("Gal", "GlcNAc", "b1-4"))
  expect_identical(inputs[[1]]$edges, list(c(1, 2), c(2, 0)))
  for (text in c(
    "Gal(a1-",
    "Gal(b1-",
    "Gal(?2-",
    "Gal(??-",
    "Gal(?1-4)GlcNAc",
    "Gal(b1-?)GlcNAc"
  )) {
    expect_error(
      sweetnet_inputs(glyparse::auto_parse(text)),
      class = "glyfp_error_unsupported"
    )
  }
  expect_identical(
    sweetnet_inputs(glyparse::auto_parse("Gal(?1-")),
    sweetnet_inputs(glyparse::auto_parse("Gal"))
  )
})

test_that("method-specific options are checked without loading torch", {
  expect_error(
    fit_fingerprint("Gal", "sweetnet", max_size = 3),
    class = "glyfp_error_argument"
  )
  expect_error(
    fit_fingerprint("Gal", "sweetnet", unknown = "ignore"),
    class = "glyfp_error_argument"
  )
  expect_error(
    fit_fingerprint("Gal", "sweetnet", unknown = "warn"),
    class = "glyfp_error_argument"
  )
  expect_false("python" %in% names(formals(fit_fingerprint)))
  expect_error(
    fit_fingerprint("Gal", checkpoint = "missing.pt"),
    class = "glyfp_error_argument"
  )
  expect_error(
    fit_fingerprint("Gal", "sweetnet", checkpoint = NULL),
    class = "glyfp_error_argument"
  )
  expect_error(
    fit_fingerprint(NA_character_, "sweetnet", na_action = "error"),
    class = "glyfp_error_missing"
  )
})

test_that("pretrained SweetNet integration preserves the public contract", {
  checkpoint <- Sys.getenv("GLYFP_TEST_CHECKPOINT")
  skip_if(
    !nzchar(checkpoint),
    "Set GLYFP_TEST_CHECKPOINT for offline pretrained integration."
  )
  skip_if_not_installed("torch")
  skip_if_not(torch::torch_is_installed(), "R torch runtime is not installed")
  x <- c(
    first = "Gal(b1-4)GlcNAc",
    branch = "Man(a1-3)[Man(a1-6)]Man",
    single = "GlcNAc",
    missing = NA_character_,
    duplicate = "Galb4GlcNAc?-"
  )
  object <- glyparse::auto_parse(x)
  e <- fit_fingerprint(x, "sweetnet", checkpoint = checkpoint)
  fp <- transform_fingerprint(x, e)
  expect_type(fp, "list")
  expect_identical(names(fp), names(x))
  expect_identical(lengths(fp), stats::setNames(rep(128L, length(x)), names(x)))
  expect_true(all(vapply(fp, is.double, TRUE)))
  expect_true(all(is.na(fp[[4]])))
  expect_identical(fp[[1]], fp[[5]])
  expect_identical(names(fp[[1]]), e$features)
  expect_identical(e$params$unknown, "error")
  expect_identical(e$params$model$device, "cpu")
  expect_identical(e$params$model$batch_size, 1L)
  expect_identical(
    e$params$model$checkpoint_sha256,
    "4ee19fc0593db0562c16585b556050694144816754750ebfe1876ed21ba0f4be"
  )
  expect_identical(
    e$params$model$vocabulary_sha256,
    "838281860e3bbddc0de91ac635bf1ba746d79bf62891e5b86ca08eff537a974e"
  )
  # Independent reference: author's SweetNet + glycans_to_emb, CPU batch_size=1.
  expect_equal(
    unname(fp[[1]][1:8]),
    c(
      -0.006824135780334473,
      0.8667784333229065,
      -1.1568934917449951,
      -0.26471179723739624,
      0.4592328667640686,
      -0.10227088630199432,
      -0.28508102893829346,
      -1.894502878189087
    ),
    tolerance = 1e-5
  )
  expect_equal(
    unname(fp[[2]][1:8]),
    c(
      -0.9553753137588501,
      0.31384938955307007,
      0.10934257507324219,
      0.46251580119132996,
      0.4480460584163666,
      -0.31087779998779297,
      -0.35662657022476196,
      -0.838372528553009
    ),
    tolerance = 1e-5
  )
  expect_identical(transform_fingerprint(object, e), fp)
  expect_identical(transform_fingerprint(x[c(5, 2, 1)], e), fp[c(5, 2, 1)])
  expect_identical(transform_fingerprint(x[1], e), fp[1])
  expect_identical(
    glycan_fingerprint(x, "sweetnet", checkpoint = checkpoint),
    fp
  )
  # Fit binds a pretrained coordinate space, not a vocabulary learned from x.
  empty <- fit_fingerprint(
    character(),
    "sweetnet",
    checkpoint = checkpoint
  )
  expect_identical(e, empty)
  expect_identical(transform_fingerprint(character(), e), list())
  expect_true(all(is.na(transform_fingerprint(NA_character_, empty)[[1]])))
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(e, path)
  expect_identical(transform_fingerprint(x[1], readRDS(path)), fp[1])
  expect_output(print(e), "CPU, 128 dimensions")
  expect_error(
    transform_fingerprint("Gal(b1-", e),
    class = "glyfp_error_unsupported"
  )
  expect_error(
    transform_fingerprint("Gal(?1-4)GlcNAc", e),
    class = "glyfp_error_unsupported"
  )
  expect_error(
    transform_fingerprint("Gal9S", e),
    class = "glyfp_error_unknown"
  )
  strict <- e
  strict$params$na_action <- "error"
  expect_error(
    transform_fingerprint(NA_character_, strict),
    class = "glyfp_error_missing"
  )
  changed <- e
  changed$params$model$readout <- "other"
  expect_error(
    transform_fingerprint("Gal", changed),
    class = "glyfp_error_backend"
  )
  changed <- e
  changed$params$model$code_sha256[["preprocessing"]] <- "changed"
  expect_error(
    transform_fingerprint("Gal", changed),
    class = "glyfp_error_backend"
  )
  changed <- e
  changed$features[1] <- "wrong"
  expect_error(
    transform_fingerprint("Gal", changed),
    class = "glyfp_error_encoder"
  )
  broken <- tempfile(fileext = ".pt")
  on.exit(unlink(broken), add = TRUE)
  writeLines("not model weights", broken)
  expect_error(
    fit_fingerprint("Gal", "sweetnet", checkpoint = broken),
    class = "glyfp_error_backend"
  )
})

test_that("native R inference agrees with the authors' model across glycan chemistry", {
  checkpoint <- Sys.getenv("GLYFP_TEST_CHECKPOINT")
  skip_if(
    !nzchar(checkpoint),
    "Set GLYFP_TEST_CHECKPOINT for pretrained parity checks"
  )
  skip_if_not_installed("torch")
  skip_if_not(torch::torch_is_installed(), "R torch runtime is not installed")
  fixture <- utils::read.csv(
    test_path("fixtures", "sweetnet-reference.csv"),
    check.names = FALSE
  )
  original_hash <- sweetnet_sha256(checkpoint)
  original_time <- file.info(checkpoint)$mtime
  threads <- torch::torch_get_num_threads()
  set.seed(725)
  rng <- .Random.seed
  encoder <- fit_fingerprint(
    fixture$glycan,
    "sweetnet",
    checkpoint = checkpoint
  )
  actual <- do.call(rbind, transform_fingerprint(fixture$glycan, encoder))
  expected <- as.matrix(fixture[-1L])
  # Absolute tolerance bounds every coordinate, including values near zero.
  expect_lt(max(abs(actual - expected)), 1e-5)
  expect_identical(sweetnet_sha256(checkpoint), original_hash)
  expect_identical(file.info(checkpoint)$mtime, original_time)
  expect_identical(torch::torch_get_num_threads(), threads)
  expect_identical(.Random.seed, rng)
  expect_identical(
    transform_fingerprint(fixture$glycan, encoder),
    transform_fingerprint(fixture$glycan, encoder)
  )
  expect_equal(
    unname(transform_fingerprint("Man(a1-6)[Man(a1-3)]Man", encoder)[[1]]),
    unname(transform_fingerprint("Man(a1-3)[Man(a1-6)]Man", encoder)[[1]]),
    tolerance = 1e-6
  )
})
