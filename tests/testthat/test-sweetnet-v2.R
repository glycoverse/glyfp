test_that("SweetNet named, local and offline paths give identical fingerprints", {
  expect_identical(formals(glycan_fingerprint)$checkpoint, "sweetnet-v2")
  for (id in c("sweetnet-v1", "sweetnet-v2")) {
    path <- Sys.getenv(
      if (id == "sweetnet-v1") {
        "GLYFP_TEST_CHECKPOINT"
      } else {
        "GLYFP_TEST_V2_CHECKPOINT"
      }
    )
    skip_if(!nzchar(path), "Set checkpoint paths for cache integration")
    cache <- tempfile()
    on.exit(unlink(cache, recursive = TRUE), add = TRUE)
    calls <- 0L
    local_mocked_bindings(sweetnet_download = function(url, destfile) {
      calls <<- calls + 1L
      stopifnot(file.copy(path, destfile))
      0L
    })
    x <- c("Gal(b1-4)GlcNAc(b1-", NA)
    local <- glycan_fingerprint(
      x,
      "sweetnet",
      checkpoint = path,
      download = FALSE
    )
    expect_identical(
      suppressMessages(glycan_fingerprint(
        x,
        "sweetnet",
        checkpoint = id,
        cache_dir = cache
      )),
      local
    )
    expect_identical(
      glycan_fingerprint(
        x,
        "sweetnet",
        checkpoint = id,
        cache_dir = cache,
        download = FALSE
      ),
      local
    )
    expect_identical(calls, 1L)
    if (id == "sweetnet-v2") {
      expect_identical(
        glycan_fingerprint(x, "sweetnet", cache_dir = cache),
        local
      )
    }
    unlink(cache, recursive = TRUE)
    expect_error(
      glycan_fingerprint(
        x,
        "sweetnet",
        checkpoint = id,
        cache_dir = cache,
        download = FALSE
      ),
      class = "glyfp_error_offline"
    )
    expect_false(dir.exists(cache))
    expect_identical(
      suppressMessages(glycan_fingerprint(
        x,
        "sweetnet",
        checkpoint = id,
        cache_dir = cache
      )),
      local
    )
    expect_identical(calls, 2L)
  }
})
