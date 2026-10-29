glycangt_test_fixture <- function() {
  skip_if_not_installed("torch")
  skip_if_not(torch::torch_is_installed(), "R torch runtime is not installed")
  skip_if_not_installed("digest")
  fixture <- readRDS(test_path("fixtures", "glycangt-untrained.rds"))
  fixture$weights <- lapply(fixture$weights, function(x) {
    torch::torch_tensor(x, dtype = torch::torch_float32(), device = "cpu")
  })
  fixture
}

test_that("GlycanGT random policies reject invalid or stale configuration", {
  for (x in list(NA, Inf, -1, 1.2, 2^32, "1", numeric())) {
    expect_error(glycangt_policy(seed = x), class = "glyfp_error_argument")
  }
  expect_error(glycangt_policy(dimensions = 0), class = "glyfp_error_argument")
  p <- glycangt_policy()
  p$algorithm <- "identity"
  expect_error(glycangt_validate_policy(p), class = "glyfp_error_backend")
  p <- glycangt_policy()
  p$batch_size <- 32L
  expect_error(glycangt_validate_policy(p), class = "glyfp_error_backend")
  # Missing checkpoints must fail before any optional backend or download.
  expect_error(
    fp_glycangt("Gal"),
    class = "glyfp_error_argument"
  )
})

test_that("native GlycanGT matches an independent untrained official model", {
  fixture <- glycangt_test_fixture()
  inputs <- glycangt_inputs(fixture$text, fixture$vocabulary)
  policy <- glycangt_policy()
  output <- glycangt_run(inputs, fixture$weights, fixture$config, policy)
  for (i in seq_along(output)) {
    expect_identical(inputs[[i]]$key, fixture$reference$key[i])
    expect_equal(
      glycangt_graph_seed(inputs[[i]]$key, policy),
      fixture$reference$seed[i]
    )
    expect_lt(
      max(abs(output[[i]] - as.double(fixture$reference[i, -(1:2)]))),
      1e-5
    )
  }
  for (n in c(1L, 3L, 8L, 65L)) {
    ref <- fixture$orf[[as.character(n)]]
    actual <- as.matrix(glycangt_orf(n, paste0("orf-case-", n), policy))
    expect_identical(dim(actual), c(n, 64L))
    expect_lt(max(abs(actual - ref)), 1e-5)
  }
})

test_that("GlycanGT results survive ordering, repeats and serialized policy", {
  fixture <- glycangt_test_fixture()
  x <- c(fixture$text[1:3], NA_character_, fixture$text[1])
  names(x) <- c("first", "branch", "sialylated", "missing", "duplicate")
  inputs <- glycangt_inputs(x, fixture$vocabulary)
  p <- glycangt_policy()
  run <- function(x, policy = p) {
    glycangt_run(x, fixture$weights, fixture$config, policy)
  }
  result <- run(inputs)
  expect_identical(names(result), names(x))
  expect_identical(result[[1]], result[[5]])
  expect_true(all(is.na(result[[4]])))
  expect_identical(run(rev(inputs)), rev(result))
  expect_identical(c(run(inputs[1:2]), run(inputs[3:5])), result)
  expect_identical(
    run(glycangt_inputs(glyparse::auto_parse(x), fixture$vocabulary)),
    result
  )
  expect_identical(
    run(glycangt_inputs(character(), fixture$vocabulary)),
    list()
  )
  file <- tempfile(fileext = ".rds")
  on.exit(unlink(file), add = TRUE)
  saveRDS(p, file)
  expect_identical(run(inputs, readRDS(file)), result)
  alternative <- glycangt_inputs("Man(a1-6)[Man(a1-3)]Man", fixture$vocabulary)
  expect_identical(run(alternative)[[1]], result[[2]])
  g <- as.list(glyparse::auto_parse(fixture$text[2]))[[1]]
  permuted <- glyrepr::glycan_structure(igraph::permute(
    g,
    rev(seq_len(igraph::vcount(g)))
  ))
  expect_identical(
    run(glycangt_inputs(permuted, fixture$vocabulary))[[1]],
    result[[2]]
  )
  expect_false(identical(
    run(inputs, glycangt_policy(seed = 1730L))[[1]],
    result[[1]]
  ))
})

test_that("GlycanGT inference preserves caller random streams", {
  fixture <- glycangt_test_fixture()
  restore <- sweetnet_restore_rng()
  on.exit(restore(), add = TRUE)
  torch::with_torch_manual_seed(
    {
      set.seed(19)
      r_seed <- .Random.seed
      torch_seed <- torch::torch_get_rng_state()$clone()
      inputs <- glycangt_inputs(fixture$text, fixture$vocabulary)
      expected <- glycangt_run(
        inputs,
        fixture$weights,
        fixture$config,
        glycangt_policy()
      )
      expect_identical(.Random.seed, r_seed)
      expect_true(torch::torch_equal(torch_seed, torch::torch_get_rng_state()))
      rm(".Random.seed", envir = globalenv())
      actual <- glycangt_run(
        inputs,
        fixture$weights,
        fixture$config,
        glycangt_policy()
      )
      expect_identical(actual, expected)
      expect_false(exists(
        ".Random.seed",
        envir = globalenv(),
        inherits = FALSE
      ))
      expect_true(torch::torch_equal(torch_seed, torch::torch_get_rng_state()))
    },
    seed = 81L
  )
})

test_that("GlycanGT trial rejects unsupported structures and wrong tensor layouts", {
  fixture <- glycangt_test_fixture()
  for (x in c("Gal", "Gal(b1-4)GlcNAc(a1-", "Gal(?1-4)GlcNAc")) {
    expect_error(
      glycangt_inputs(x, fixture$vocabulary),
      class = "glyfp_error_unsupported"
    )
  }
  expect_error(
    glycangt_inputs("Gal(b1-3)GlcNAc", fixture$vocabulary),
    class = "glyfp_error_unknown"
  )
  inputs <- glycangt_inputs(fixture$text[1], fixture$vocabulary)
  bad <- fixture$weights
  bad[["encoder.graph_feature.orf_encoder.weight"]] <- NULL
  expect_error(
    glycangt_run(inputs, bad, fixture$config, glycangt_policy()),
    class = "glyfp_error_backend"
  )
  bad <- fixture$config
  bad$heads <- 3L
  expect_error(
    glycangt_run(inputs, fixture$weights, bad, glycangt_policy()),
    class = "glyfp_error_backend"
  )
  expect_error(
    glycangt_run(
      inputs,
      fixture$weights,
      fixture$config,
      glycangt_policy(dimensions = 32L)
    ),
    class = "glyfp_error_backend"
  )
})


test_that("GlycanGT fixes float32 independent of the caller default", {
  fixture <- glycangt_test_fixture()
  inputs <- glycangt_inputs(fixture$text[1], fixture$vocabulary)
  expected <- glycangt_run(
    inputs,
    fixture$weights,
    fixture$config,
    glycangt_policy()
  )
  old <- torch::torch_get_default_dtype()
  on.exit(torch::torch_set_default_dtype(old), add = TRUE)
  torch::torch_set_default_dtype(torch::torch_float64())
  expect_identical(
    glycangt_run(inputs, fixture$weights, fixture$config, glycangt_policy()),
    expected
  )
  bad <- fixture$config
  bad$heads <- NA_integer_
  expect_error(
    glycangt_run(inputs, fixture$weights, bad, glycangt_policy()),
    class = "glyfp_error_backend"
  )
})
