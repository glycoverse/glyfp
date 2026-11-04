test_that("GIFFLAR defaults and explicit cache options reach its backend", {
  calls <- list()
  local_mocked_bindings(gifflar_backend = function(
    glycans,
    checkpoint,
    cache_dir,
    download
  ) {
    calls[[length(calls) + 1L]] <<- list(checkpoint, cache_dir, download)
    lapply(glycans, function(x) numeric(384))
  })
  fp_gifflar("GlcNAc(b1-")
  fp_gifflar("GlcNAc(b1-", "local.pt", cache_dir = "custom", download = FALSE)
  expect_identical(
    calls[[1]],
    list("gifflar-species", tools::R_user_dir("glyfp", "cache"), TRUE)
  )
  expect_identical(calls[[2]], list("local.pt", "custom", FALSE))
})

test_that("GIFFLAR uses pinned verified downloads and offline caches", {
  payload <- charToRaw("fake GIFFLAR weights")
  hash <- digest::digest(payload, algo = "sha256", serialize = FALSE)
  spec <- gifflar_spec()
  spec$sha256 <- hash
  cache <- tempfile()
  on.exit(unlink(cache, recursive = TRUE), add = TRUE)
  calls <- 0L
  local_mocked_bindings(
    gifflar_spec = function() spec,
    sweetnet_download = function(url, destfile) {
      calls <<- calls + 1L
      expect_identical(
        url,
        paste0(
          "https://huggingface.co/DBojar/glycowork_models/resolve/",
          spec$checkpoint_revision,
          "/glycowork_gifflar_species.pt"
        )
      )
      writeBin(payload, destfile)
      0L
    }
  )
  expect_error(
    gifflar_resolve("gifflar-species", cache, FALSE),
    class = "glyfp_error_offline"
  )
  expect_false(dir.exists(cache))
  first <- suppressMessages(gifflar_resolve("gifflar-species", cache, TRUE))
  expect_identical(first$id, "gifflar-species")
  expect_identical(gifflar_resolve("gifflar-species", cache, FALSE), first)
  expect_identical(gifflar_resolve("gifflar-species", cache, TRUE), first)
  expect_identical(calls, 1L)
  expect_identical(gifflar_resolve(first$path, NULL, FALSE)$path, first$path)
  writeLines("corrupt", first$path)
  expect_error(
    gifflar_resolve("gifflar-species", cache, TRUE),
    class = "glyfp_error_checksum"
  )
  expect_error(
    gifflar_resolve(first$path, cache, TRUE),
    class = "glyfp_error_checksum"
  )
  expect_identical(calls, 1L)
  expect_length(list.files(cache, pattern = "[.]part$"), 0L)
})

test_that("GIFFLAR validates options and cleans failed downloads", {
  cache <- tempfile()
  on.exit(unlink(cache, recursive = TRUE), add = TRUE)
  for (bad in list(
    NULL,
    NA_character_,
    character(),
    "",
    "https://example.org/model.pt"
  )) {
    expect_error(
      gifflar_resolve(bad, cache, TRUE),
      class = "glyfp_error_argument"
    )
  }
  expect_error(
    gifflar_resolve("gifflar-species", cache, NA),
    class = "glyfp_error_argument"
  )
  expect_error(
    gifflar_resolve("gifflar-species", NULL, TRUE),
    class = "glyfp_error_argument"
  )
  for (failure in c("error", "status", "checksum")) {
    local_mocked_bindings(sweetnet_download = function(url, destfile) {
      writeLines("partial", destfile)
      if (failure == "error") {
        stop("connection lost")
      }
      if (failure == "status") {
        return(1L)
      }
      0L
    })
    expect_error(
      suppressMessages(gifflar_resolve("gifflar-species", cache, TRUE)),
      class = if (failure == "checksum") {
        "glyfp_error_checksum"
      } else {
        "glyfp_error_download"
      }
    )
    expect_length(list.files(cache, all.files = TRUE, no.. = TRUE), 0L)
  }
})
