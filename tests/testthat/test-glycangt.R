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
        y <- glyfp:::glycangt_backend(x, checkpoint)
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
