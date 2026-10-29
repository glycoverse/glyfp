test_that("GIFFLAR requires an explicit local checkpoint without network", {
  local_mocked_bindings(
    download.file = function(...) stop("Unexpected network"),
    .package = "utils"
  )
  expect_error(
    fit_fingerprint("Gal", "gifflar"),
    class = "glyfp_error_argument"
  )
  for (path in list(
    NULL,
    character(),
    NA_character_,
    "",
    "https://example.org/model.pt"
  )) {
    expect_error(
      fit_fingerprint("Gal", "gifflar", checkpoint = path),
      class = "glyfp_error_argument"
    )
  }
  expect_error(
    fit_fingerprint("Gal", "gifflar", checkpoint = tempfile()),
    class = "glyfp_error_backend"
  )
  for (option in list(
    list(download = TRUE),
    list(download = NA),
    list(cache_dir = tempfile()),
    list(max_size = 2),
    list(unknown = "ignore")
  )) {
    expect_error(
      do.call(
        fit_fingerprint,
        c(list("Gal", "gifflar", checkpoint = "local.pt"), option)
      ),
      class = "glyfp_error_argument"
    )
  }
})

test_that("native GIFFLAR molecular graphs match all fixed residue templates", {
  skip_if_not_installed("digest")
  templates <- gifflar_templates()
  vocabulary <- sweetnet_vocabulary()
  references <- readRDS(test_path("fixtures", "gifflar-reference.rds"))
  sorted <- function(x) {
    if (!nrow(x)) {
      return(x)
    }
    unname(x[do.call(order, as.data.frame(x)), , drop = FALSE])
  }
  for (reference in references) {
    object <- glyrepr::as_glycan_structure(reference$glycan)
    fragment <- gifflar_smiles(as.list(object)[[1]], templates)
    expect_identical(fragment$text, reference$smiles, info = reference$glycan)
    input <- gifflar_inputs(object, templates, vocabulary)[[1]]
    expect_equal(input$x, reference$input$x, info = reference$glycan)
    for (relation in names(reference$input$edges)) {
      actual <- input$edges[[relation]]
      expected <- reference$input$edges[[relation]]
      if (!is.null(reference$input$edge_attr[[relation]])) {
        actual <- cbind(actual, input$edge_attr[[relation]])
        expected <- cbind(expected, reference$input$edge_attr[[relation]])
      }
      expect_equal(
        sorted(actual),
        sorted(expected),
        info = paste(reference$glycan, relation)
      )
    }
    expect_identical(
      gifflar_inputs(reference$canonical, templates, vocabulary)[[1]],
      input
    )
  }
  expect_identical(gifflar_inputs(character(), templates, vocabulary), list())
  expect_identical(
    gifflar_inputs(NA_character_, templates, vocabulary),
    list(NULL)
  )
  for (bad in c(
    "Galf",
    "Gal(?1-4)GlcNAc",
    "Gal(b1-4)GlcNAc(b1-",
    "Gal(b1-1)Glc",
    "Gal(b2-4)Glc",
    "Gal(b1-2)GlcN",
    "Gal(b1-6)Glc6P",
    "Neu5Ac(a1-3)Gal"
  )) {
    expect_error(
      gifflar_inputs(glyrepr::as_glycan_structure(bad), templates, vocabulary),
      class = "glyfp_error_unsupported",
      info = bad
    )
  }
})

test_that("native GIFFLAR network matches an independent small untrained model", {
  skip_if_not_installed("torch")
  skip_if_not_installed("digest")
  skip_if_not(torch::torch_is_installed(), "R torch runtime is not installed")
  toy <- readRDS(test_path("fixtures", "gifflar-untrained.rds"))
  weights <- lapply(toy$weights, function(x) {
    torch::torch_tensor(x, dtype = torch::torch_float32(), device = "cpu")
  })
  for (reference in toy$reference) {
    input <- gifflar_inputs(
      reference$glycan,
      gifflar_templates(),
      sweetnet_vocabulary()
    )[[1]]
    expect_equal(
      gifflar_infer(input, weights, layers = 2L),
      reference$embedding,
      tolerance = 1e-5
    )
  }
})

test_that("pretrained native GIFFLAR preserves all encoder and batch contracts", {
  checkpoint <- Sys.getenv("GLYFP_TEST_GIFFLAR_CHECKPOINT")
  skip_if(
    !nzchar(checkpoint),
    "Set GLYFP_TEST_GIFFLAR_CHECKPOINT for real integration"
  )
  skip_if_not_installed("torch")
  skip_if_not(torch::torch_is_installed(), "R torch runtime is not installed")
  local_mocked_bindings(
    download.file = function(...) stop("Unexpected network"),
    .package = "utils"
  )
  references <- readRDS(test_path("fixtures", "gifflar-reference.rds"))
  glycans <- vapply(references, `[[`, "", "canonical")
  before <- sweetnet_sha256(checkpoint)
  encoder <- fit_fingerprint(glycans, "gifflar", checkpoint = checkpoint)
  expect_identical(encoder$params$model$spec$sha256, before)
  expect_identical(
    encoder$params$model$spec$vocabulary_sha256,
    sweetnet_sha256(system.file(
      "extdata",
      "sweetnet-vocabulary.tsv",
      package = "glyfp"
    ))
  )
  actual <- transform_fingerprint(glycans, encoder)
  expected <- do.call(rbind, lapply(references, `[[`, "embedding"))
  expect_lt(max(abs(do.call(rbind, actual) - expected)), 1e-5)
  x <- c(
    a = glycans[1],
    b = glycans[2],
    missing = NA_character_,
    duplicate = glycans[1]
  )
  output <- transform_fingerprint(x, encoder)
  expect_identical(names(output), names(x))
  expect_identical(unname(lengths(output)), rep(384L, 4L))
  expect_identical(output[[1]], output[[4]])
  expect_true(all(is.na(output[[3]])))
  expect_type(output[[1]], "double")
  expect_identical(names(output[[1]]), encoder$features)
  expect_identical(transform_fingerprint(rev(x), encoder), rev(output))
  expect_identical(
    c(
      transform_fingerprint(x[1:2], encoder),
      transform_fingerprint(x[3:4], encoder)
    ),
    output
  )
  expect_identical(
    transform_fingerprint(glyparse::auto_parse(x), encoder),
    output
  )
  expect_identical(
    glycan_fingerprint(x, "gifflar", checkpoint = checkpoint, download = FALSE),
    output
  )
  # Canonical branch storage must not affect the molecular model.
  branched <- c("Man(a1-3)[Man(a1-6)]Man", "Man(a1-6)[Man(a1-3)]Man")
  branched <- transform_fingerprint(branched, encoder)
  expect_identical(branched[[1]], branched[[2]])
  saved <- tempfile(fileext = ".rds")
  corrupt <- tempfile(fileext = ".pt")
  on.exit(unlink(c(saved, corrupt)), add = TRUE)
  saveRDS(encoder, saved)
  expect_lt(file.info(saved)$size, 20000)
  expect_identical(transform_fingerprint(x, readRDS(saved)), output)
  empty <- fit_fingerprint(character(), "gifflar", checkpoint = checkpoint)
  expect_identical(empty$features, encoder$features)
  expect_identical(transform_fingerprint(character(), empty), list())
  missing <- fit_fingerprint(NA_character_, "gifflar", checkpoint = checkpoint)
  expect_true(all(is.na(transform_fingerprint(NA_character_, missing)[[1]])))
  expect_error(
    fit_fingerprint(
      NA_character_,
      "gifflar",
      checkpoint = checkpoint,
      na_action = "error"
    ),
    class = "glyfp_error_missing"
  )
  strict <- fit_fingerprint(
    glycans[1],
    "gifflar",
    checkpoint = checkpoint,
    na_action = "error"
  )
  expect_error(
    transform_fingerprint(NA_character_, strict),
    class = "glyfp_error_missing"
  )
  changed <- encoder
  changed$params$model$batch_size <- 2L
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
  writeLines("corrupted", corrupt)
  expect_error(
    fit_fingerprint(glycans[1], "gifflar", checkpoint = corrupt),
    class = "glyfp_error_backend"
  )
  changed <- encoder
  changed$params$checkpoint <- corrupt
  expect_error(
    transform_fingerprint(character(), changed),
    class = "glyfp_error_backend"
  )
  expect_identical(sweetnet_sha256(checkpoint), before)
  expect_output(print(encoder), "384 dimensions")
})

test_that("first GIFFLAR load preserves RNG and never installs a runtime", {
  checkpoint <- Sys.getenv("GLYFP_TEST_GIFFLAR_CHECKPOINT")
  skip_if(
    !nzchar(checkpoint),
    "Set GLYFP_TEST_GIFFLAR_CHECKPOINT for first-load integration"
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
          set.seed(684)
          before <- .Random.seed
          Sys.setenv(TORCH_INSTALL = "1")
        } else {
          if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
            rm(".Random.seed", envir = globalenv())
          }
          Sys.unsetenv("TORCH_INSTALL")
        }
        result <- glyfp:::gifflar_backend(x, checkpoint)$embeddings
        list(
          output = length(result[[1]]) == 384L && all(is.finite(result[[1]])),
          rng = if (seeded) {
            identical(before, .Random.seed)
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
