test_that("GlycanGT requires an explicit local file and cannot download", {
  local_mocked_bindings(
    download.file = function(...) stop("Unexpected network call"),
    .package = "utils"
  )
  cache <- tempfile("unused-glycangt-cache-")
  expect_error(
    fit_fingerprint("Gal(b1-4)GlcNAc", "glycangt"),
    class = "glyfp_error_argument"
  )
  for (file in list(
    NULL,
    NA_character_,
    character(),
    "",
    "https://huggingface.co/model.pt"
  )) {
    expect_error(
      fit_fingerprint("Gal(b1-4)GlcNAc", "glycangt", checkpoint = file),
      class = "glyfp_error_argument"
    )
  }
  expect_error(
    fit_fingerprint("Gal(b1-4)GlcNAc", "glycangt", checkpoint = tempfile()),
    class = "glyfp_error_backend"
  )
  expect_error(
    fit_fingerprint(
      "Gal(b1-4)GlcNAc",
      "glycangt",
      checkpoint = "local.pt",
      download = TRUE
    ),
    class = "glyfp_error_argument"
  )
  expect_error(
    fit_fingerprint(
      "Gal(b1-4)GlcNAc",
      "glycangt",
      checkpoint = "local.pt",
      cache_dir = cache
    ),
    class = "glyfp_error_argument"
  )
  expect_false(dir.exists(cache))
  expect_error(
    fit_fingerprint(
      "Gal(b1-4)GlcNAc",
      "glycangt",
      checkpoint = "local.pt",
      unknown = "ignore"
    ),
    class = "glyfp_error_argument"
  )
  expect_error(
    fit_fingerprint(
      "Gal(b1-4)GlcNAc",
      "glycangt",
      checkpoint = "local.pt",
      max_size = 2
    ),
    class = "glyfp_error_argument"
  )
})

test_that("real GlycanGT inference matches the official large model", {
  path <- Sys.getenv("GLYFP_TEST_GLYCANGT_CHECKPOINT")
  skip_if(
    !nzchar(path),
    "Set GLYFP_TEST_GLYCANGT_CHECKPOINT for local pretrained integration"
  )
  skip_if_not_installed("torch")
  skip_if_not(torch::torch_is_installed(), "R torch runtime is not installed")
  local_mocked_bindings(
    download.file = function(...) stop("Unexpected network call"),
    .package = "utils"
  )
  reference <- read.csv(
    test_path("fixtures", "glycangt-reference.csv"),
    check.names = FALSE
  )
  before <- sweetnet_sha256(path)
  encoder <- fit_fingerprint(reference$key, "glycangt", checkpoint = path)
  expect_identical(encoder$params$model$checkpoint_sha256, before)
  expect_identical(encoder$params$model$config$dimensions, 768L)
  expect_identical(encoder$params$model$policy, glycangt_policy())
  actual <- transform_fingerprint(reference$key, encoder)
  expect_lt(
    max(abs(do.call(rbind, actual) - as.matrix(reference[, -(1:2)]))),
    1e-5
  )
  x <- c(
    a = reference$key[1],
    b = reference$key[2],
    missing = NA_character_,
    duplicate = reference$key[1]
  )
  output <- transform_fingerprint(x, encoder)
  expect_type(output, "list")
  expect_identical(names(output), names(x))
  expect_identical(unname(lengths(output)), rep(768L, 4L))
  expect_identical(output[[1]], output[[4]])
  expect_true(all(is.na(output[[3]])))
  expect_identical(names(output[[1]]), encoder$features)
  expect_identical(
    transform_fingerprint(glyparse::auto_parse(x), encoder),
    output
  )
  expect_identical(transform_fingerprint(rev(x), encoder), rev(output))
  expect_identical(
    c(
      transform_fingerprint(x[1:2], encoder),
      transform_fingerprint(x[3:4], encoder)
    ),
    output
  )
  expect_identical(
    glycan_fingerprint(x, "glycangt", checkpoint = path, download = FALSE),
    output
  )
  saved <- tempfile(fileext = ".rds")
  corrupt <- tempfile(fileext = ".pt")
  on.exit(unlink(c(saved, corrupt)), add = TRUE)
  saveRDS(encoder, saved)
  expect_lt(file.info(saved)$size, 20000)
  expect_identical(transform_fingerprint(x, readRDS(saved)), output)
  empty <- fit_fingerprint(character(), "glycangt", checkpoint = path)
  expect_identical(empty$features, encoder$features)
  expect_identical(transform_fingerprint(character(), empty), list())
  missing <- fit_fingerprint(
    c(NA_character_, NA_character_),
    "glycangt",
    checkpoint = path
  )
  expect_identical(
    transform_fingerprint(NA_character_, missing)[[1]],
    stats::setNames(rep(NA_real_, 768L), missing$features)
  )
  expect_error(
    fit_fingerprint(
      NA_character_,
      "glycangt",
      checkpoint = path,
      na_action = "error"
    ),
    class = "glyfp_error_missing"
  )
  strict <- fit_fingerprint(
    reference$key[1],
    "glycangt",
    checkpoint = path,
    na_action = "error"
  )
  expect_error(
    transform_fingerprint(NA_character_, strict),
    class = "glyfp_error_missing"
  )
  changed <- encoder
  changed$params$model$policy$seed <- 1L
  expect_error(
    transform_fingerprint(character(), changed),
    class = "glyfp_error_backend"
  )
  changed <- encoder
  changed$features <- changed$features[-1]
  expect_error(
    transform_fingerprint(character(), changed),
    class = "glyfp_error_encoder"
  )
  writeLines("not a checkpoint", corrupt)
  expect_error(
    fit_fingerprint(reference$key[1], "glycangt", checkpoint = corrupt),
    class = "glyfp_error_backend"
  )
  changed <- encoder
  changed$params$checkpoint <- corrupt
  # A warm tensor cache must not skip hashing the supplied local artifact.
  expect_error(
    transform_fingerprint(character(), changed),
    class = "glyfp_error_backend"
  )
  expect_identical(sweetnet_sha256(path), before)
  expect_output(print(encoder), "768 dimensions")
})

test_that("first GlycanGT backend load preserves RNG and installation policy", {
  checkpoint <- Sys.getenv("GLYFP_TEST_GLYCANGT_CHECKPOINT")
  skip_if(
    !nzchar(checkpoint),
    "Set GLYFP_TEST_GLYCANGT_CHECKPOINT for first-load integration"
  )
  skip_if_not_installed("torch")
  for (seeded in c(TRUE, FALSE)) {
    result <- callr::r(
      function(package, checkpoint, seeded) {
        if (file.exists(file.path(package, "R", "fingerprint.R"))) {
          pkgload::load_all(package, quiet = TRUE)
        } else {
          library("glyfp", lib.loc = dirname(package), character.only = TRUE)
        }
        stopifnot(!"torch" %in% loadedNamespaces())
        x <- glyparse::auto_parse("Gal(b1-4)GlcNAc")
        if (seeded) {
          set.seed(725)
          before <- .Random.seed
          Sys.setenv(TORCH_INSTALL = "1")
        } else {
          if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
            rm(".Random.seed", envir = globalenv())
          }
          Sys.unsetenv("TORCH_INSTALL")
        }
        y <- glyfp:::glycangt_backend(x, checkpoint)$embeddings
        list(
          output = length(y) == 1L &&
            length(y[[1]]) == 768L &&
            all(is.finite(y[[1]])),
          rng = if (seeded) {
            identical(.Random.seed, before)
          } else {
            !exists(".Random.seed", envir = globalenv(), inherits = FALSE)
          },
          install = if (seeded) {
            identical(Sys.getenv("TORCH_INSTALL"), "1")
          } else {
            is.na(Sys.getenv("TORCH_INSTALL", unset = NA_character_))
          }
        )
      },
      args = list(
        package = getNamespaceInfo("glyfp", "path"),
        checkpoint = checkpoint,
        seeded = seeded
      ),
      libpath = .libPaths()
    )
    expect_identical(result, list(output = TRUE, rng = TRUE, install = TRUE))
  }
})
