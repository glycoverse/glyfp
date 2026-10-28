test_that("SweetNet version identity and download policy are validated", {
  expect_identical(formals(fit_fingerprint)$checkpoint, "sweetnet-v2")
  expect_error(
    fit_fingerprint("Gal", "sweetnet", download = NA),
    class = "glyfp_error_argument"
  )
  expect_error(
    fit_fingerprint("Gal", "sweetnet", cache_dir = NULL),
    class = "glyfp_error_argument"
  )
  expect_error(
    fit_fingerprint(
      "Gal",
      "sweetnet",
      checkpoint = c("sweetnet-v1", "sweetnet-v2")
    ),
    class = "glyfp_error_argument"
  )
  catalog <- sweetnet_checkpoints()
  expect_identical(names(catalog), c("sweetnet-v1", "sweetnet-v2"))
  expect_false(identical(catalog[[1]]$sha256, catalog[[2]]$sha256))
})

test_that("the same legacy model is reusable by name, local path, and offline encoder", {
  checkpoint <- Sys.getenv("GLYFP_TEST_CHECKPOINT")
  skip_if(
    !nzchar(checkpoint),
    "Set GLYFP_TEST_CHECKPOINT for pretrained cache integration"
  )
  skip_if_not_installed("torch")
  skip_if_not(torch::torch_is_installed(), "R torch runtime is not installed")
  cache <- tempfile()
  on.exit(unlink(cache, recursive = TRUE), add = TRUE)
  calls <- 0L
  local_mocked_bindings(sweetnet_download = function(url, destfile) {
    calls <<- calls + 1L
    stopifnot(file.copy(checkpoint, destfile))
    0L
  })
  a <- c("Gal(b1-4)GlcNAc", "GlcNAc", NA_character_)
  named <- suppressMessages(fit_fingerprint(
    a,
    "sweetnet",
    checkpoint = "sweetnet-v1",
    cache_dir = cache
  ))
  offline <- fit_fingerprint(
    a,
    "sweetnet",
    checkpoint = "sweetnet-v1",
    cache_dir = cache,
    download = FALSE
  )
  local <- fit_fingerprint(
    a,
    "sweetnet",
    checkpoint = checkpoint,
    download = FALSE
  )
  expect_identical(named$params$model$model, "sweetnet-v1")
  expect_identical(named$params$model, local$params$model)
  expect_identical(
    transform_fingerprint(a, named),
    transform_fingerprint(a, offline)
  )
  expect_identical(
    transform_fingerprint(a, local),
    transform_fingerprint(a, offline)
  )
  expect_identical(calls, 1L)
  rds <- tempfile()
  on.exit(unlink(rds), add = TRUE)
  saveRDS(offline, rds)
  expect_identical(
    transform_fingerprint(a, readRDS(rds)),
    transform_fingerprint(a, local)
  )
  unlink(cache, recursive = TRUE)
  stale <- named
  stale$params$model$readout <- "changed"
  expect_error(transform_fingerprint(a, stale), class = "glyfp_error_backend")
  expect_identical(calls, 1L)
  expect_error(transform_fingerprint(a, offline), class = "glyfp_error_offline")
  expect_identical(
    suppressMessages(transform_fingerprint(a, named)),
    transform_fingerprint(a, local)
  )
  expect_identical(calls, 2L)
  changed <- named
  changed$params$checkpoint <- "sweetnet-v2"
  expect_error(transform_fingerprint(a, changed), class = "glyfp_error_encoder")
})

test_that("SweetNet v2 matches the official pretrained model", {
  checkpoint <- Sys.getenv("GLYFP_TEST_V2_CHECKPOINT")
  skip_if(
    !nzchar(checkpoint),
    "Set GLYFP_TEST_V2_CHECKPOINT for actual v2 parity"
  )
  skip_if_not_installed("torch")
  skip_if_not(torch::torch_is_installed(), "R torch runtime is not installed")
  fixture <- utils::read.csv(
    test_path("fixtures", "sweetnet-v2-reference.csv"),
    check.names = FALSE
  )
  original <- sweetnet_sha256(checkpoint)
  encoder <- fit_fingerprint(
    fixture$glycan,
    "sweetnet",
    checkpoint = checkpoint,
    download = FALSE
  )
  fp <- transform_fingerprint(fixture$glycan, encoder)
  expect_identical(encoder$params$model$model, "sweetnet-v2")
  expect_identical(
    encoder$params$model$checkpoint_sha256,
    sweetnet_checkpoints()[["sweetnet-v2"]]$sha256
  )
  expect_lt(max(abs(do.call(rbind, fp) - as.matrix(fixture[-1L]))), 1e-5)
  expect_identical(sweetnet_sha256(checkpoint), original)
  expect_identical(
    transform_fingerprint(glyparse::auto_parse(fixture$glycan), encoder),
    fp
  )
  expect_identical(
    transform_fingerprint(fixture$glycan[c(3, 1, 3)], encoder),
    fp[c(3, 1, 3)]
  )
  expect_identical(transform_fingerprint(character(), encoder), list())
  expect_true(all(is.na(transform_fingerprint(NA_character_, encoder)[[1L]])))
  cache <- tempfile()
  on.exit(unlink(cache, recursive = TRUE), add = TRUE)
  local_mocked_bindings(sweetnet_download = function(url, destfile) {
    stopifnot(file.copy(checkpoint, destfile))
    0L
  })
  default <- suppressMessages(fit_fingerprint(
    fixture$glycan,
    "sweetnet",
    cache_dir = cache
  ))
  expect_identical(default$params$model, encoder$params$model)
  expect_identical(transform_fingerprint(fixture$glycan, default), fp)
})
