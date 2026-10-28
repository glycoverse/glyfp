fake_checkpoint_cache <- function(env = parent.frame()) {
  cache <- tempfile("glyfp-test-cache-")
  payload <- charToRaw("test-only checkpoint bytes\n")
  hash <- digest::digest(payload, algo = "sha256", serialize = FALSE)
  spec <- list(revision = "fixed-revision", file = "weights.pt", sha256 = hash)
  local_mocked_bindings(
    sweetnet_checkpoints = function() {
      list(`sweetnet-v1` = spec, `sweetnet-v2` = spec)
    },
    sweetnet_download = function(...) {
      stop("Network must be mocked in cache tests")
    },
    .env = env
  )
  list(dir = cache, payload = payload, hash = hash)
}

test_that("named checkpoints download once and publish only verified files", {
  fake <- fake_checkpoint_cache()
  on.exit(unlink(fake$dir, recursive = TRUE), add = TRUE)
  calls <- 0L
  local_mocked_bindings(sweetnet_download = function(url, destfile) {
    calls <<- calls + 1L
    expect_identical(
      url,
      "https://huggingface.co/DBojar/glycowork_models/resolve/fixed-revision/weights.pt"
    )
    expect_identical(dirname(destfile), fake$dir)
    expect_false(file.exists(file.path(
      fake$dir,
      paste0("sweetnet-v2-", fake$hash, ".pt")
    )))
    writeBin(fake$payload, destfile)
    0L
  })
  first <- suppressMessages(sweetnet_resolve("sweetnet-v2", fake$dir, TRUE))
  expect_identical(first$id, "sweetnet-v2")
  expect_identical(first$source, "sweetnet-v2")
  expect_identical(sweetnet_sha256(first$path), fake$hash)
  expect_identical(sweetnet_resolve("sweetnet-v2", fake$dir, FALSE), first)
  expect_identical(sweetnet_resolve("sweetnet-v2", fake$dir, TRUE), first)
  expect_identical(calls, 1L)
  expect_length(list.files(fake$dir, pattern = "\\.part$"), 0L)
})

test_that("offline misses and damaged caches never download", {
  fake <- fake_checkpoint_cache()
  on.exit(unlink(fake$dir, recursive = TRUE), add = TRUE)
  expect_error(
    sweetnet_resolve("sweetnet-v2", fake$dir, FALSE),
    class = "glyfp_error_offline"
  )
  expect_false(dir.exists(fake$dir))
  dir.create(fake$dir)
  path <- file.path(fake$dir, paste0("sweetnet-v2-", fake$hash, ".pt"))
  writeLines("damaged", path)
  for (allow in c(TRUE, FALSE)) {
    expect_error(
      sweetnet_resolve("sweetnet-v2", fake$dir, allow),
      class = "glyfp_error_checksum"
    )
  }
  expect_identical(readLines(path), "damaged")
})

test_that("failed or corrupt downloads leave no final or partial file", {
  for (failure in c("error", "status", "checksum")) {
    fake <- fake_checkpoint_cache()
    local_mocked_bindings(sweetnet_download = function(url, destfile) {
      writeLines("partial response", destfile)
      if (failure == "error") {
        stop("connection lost")
      }
      if (failure == "status") {
        return(1L)
      }
      0L
    })
    expect_error(
      suppressMessages(sweetnet_resolve("sweetnet-v2", fake$dir, TRUE)),
      class = if (failure == "checksum") {
        "glyfp_error_checksum"
      } else {
        "glyfp_error_download"
      }
    )
    expect_length(list.files(fake$dir, all.files = TRUE, no.. = TRUE), 0L)
    unlink(fake$dir, recursive = TRUE)
  }
})

test_that("explicit local files take precedence over caching and are validated", {
  fake <- fake_checkpoint_cache()
  path <- tempfile(fileext = ".pt")
  on.exit(unlink(c(path, fake$dir), recursive = TRUE), add = TRUE)
  # Give only one registered model this hash so local identification is unambiguous.
  specs <- sweetnet_checkpoints()
  specs[["sweetnet-v1"]]$sha256 <- "different"
  local_mocked_bindings(sweetnet_checkpoints = function() specs)
  writeBin(fake$payload, path)
  result <- sweetnet_resolve(path, fake$dir, FALSE)
  expect_identical(result$id, "sweetnet-v2")
  expect_identical(result$path, normalizePath(path))
  expect_identical(result$source, normalizePath(path))
  expect_false(dir.exists(fake$dir))
  writeLines("corrupt", path)
  expect_error(
    sweetnet_resolve(path, fake$dir, TRUE),
    class = "glyfp_error_checksum"
  )
  expect_false(dir.exists(fake$dir))
})

test_that("a concurrent verified cache winner is reused", {
  fake <- fake_checkpoint_cache()
  on.exit(unlink(fake$dir, recursive = TRUE), add = TRUE)
  local_mocked_bindings(sweetnet_download = function(url, destfile) {
    writeBin(fake$payload, destfile)
    writeBin(
      fake$payload,
      file.path(fake$dir, paste0("sweetnet-v2-", fake$hash, ".pt"))
    )
    0L
  })
  result <- suppressMessages(sweetnet_resolve("sweetnet-v2", fake$dir, TRUE))
  expect_identical(sweetnet_sha256(result$path), fake$hash)
  expect_length(list.files(fake$dir), 1L)
})

test_that("checkpoint options and cache failures are explicit", {
  fake <- fake_checkpoint_cache()
  on.exit(unlink(fake$dir, recursive = TRUE), add = TRUE)
  for (bad in list(NULL, NA_character_, character(), c("a", "b"))) {
    expect_error(
      sweetnet_resolve(bad, fake$dir, TRUE),
      class = "glyfp_error_argument"
    )
  }
  expect_error(
    sweetnet_resolve("sweetnet-v2", fake$dir, NA),
    class = "glyfp_error_argument"
  )
  expect_error(sweetnet_cache_dir(NULL), class = "glyfp_error_argument")
  writeLines("file blocks directory", fake$dir)
  expect_error(
    sweetnet_resolve("sweetnet-v2", fake$dir, TRUE),
    class = "glyfp_error_cache"
  )
})

test_that("concurrent first directory creation is accepted", {
  fake <- fake_checkpoint_cache()
  on.exit(unlink(fake$dir, recursive = TRUE), add = TRUE)
  local_mocked_bindings(
    sweetnet_create_cache = function(path) {
      dir.create(path)
      FALSE # another process created it between the existence check and mkdir
    },
    sweetnet_download = function(url, destfile) {
      writeBin(fake$payload, destfile)
      0L
    }
  )
  result <- suppressMessages(sweetnet_resolve("sweetnet-v2", fake$dir, TRUE))
  expect_identical(sweetnet_sha256(result$path), fake$hash)
})
