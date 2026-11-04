test_that("optional dependency failures preserve RNG and suppress automatic installation", {
  install <- Sys.getenv("TORCH_INSTALL", unset = NA_character_)
  on.exit(
    {
      if (is.na(install)) {
        Sys.unsetenv("TORCH_INSTALL")
      } else {
        Sys.setenv(TORCH_INSTALL = install)
      }
    },
    add = TRUE
  )
  Sys.setenv(TORCH_INSTALL = "1")
  set.seed(725)
  before <- .Random.seed
  local_mocked_bindings(sweetnet_require_namespace = function(pkg) {
    expect_identical(Sys.getenv("TORCH_INSTALL"), "0")
    # Model the RNG consumption in torch's namespace startup before a failure.
    sample.int(100000L, 1L)
    FALSE
  })
  expect_error(sweetnet_dependencies(), class = "glyfp_error_backend")
  expect_identical(Sys.getenv("TORCH_INSTALL"), "1")
  expect_identical(.Random.seed, before)
})

test_that("structural methods do not load torch in a fresh R process", {
  result <- callr::r(
    function(package) {
      if (file.exists(file.path(package, "R", "fingerprint.R"))) {
        pkgload::load_all(package, quiet = TRUE)
      } else {
        library("glyfp", lib.loc = dirname(package), character.only = TRUE)
      }
      stopifnot(!"torch" %in% loadedNamespaces())
      # Isolate optional backend startup from igraph, which initializes an absent
      # R RNG state during graph parsing. Build the graph before clearing it.
      inputs <- glyfp:::sweetnet_inputs(glyparse::auto_parse("Gal"))
      glyfp::fp_motif_binary("Gal(b1-")
      glyfp::fp_motif_count("Gal(b1-")
      glyfp::fp_circular_binary("Gal(b1-")
      glyfp::fp_residue_pair("Gal(b1-")
      !"torch" %in% loadedNamespaces()
    },
    args = list(package = getNamespaceInfo("glyfp", "path")),
    libpath = .libPaths()
  )
  expect_true(result)
})

test_that("first native inference preserves present and absent caller state", {
  checkpoint <- Sys.getenv("GLYFP_TEST_CHECKPOINT")
  skip_if(
    !nzchar(checkpoint),
    "Set GLYFP_TEST_CHECKPOINT for first-load integration"
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
        # Isolate optional backend startup from igraph, which initializes an absent
        # R RNG state during graph parsing. Build the graph before clearing it.
        inputs <- glyfp:::sweetnet_inputs(glyparse::auto_parse("Gal"))
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
        fingerprint <- glyfp:::sweetnet_run(inputs, checkpoint)
        list(
          output = length(fingerprint) == 1L &&
            length(fingerprint[[1L]]) == 128L &&
            all(is.finite(fingerprint[[1L]])),
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
