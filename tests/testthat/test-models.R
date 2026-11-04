test_that("model options and local-only checkpoints are enforced", {
  local_mocked_bindings(sweetnet_dependencies = function() NULL)
  local_mocked_bindings(
    download.file = function(...) stop("Unexpected network"),
    .package = "utils"
  )
  x <- "Gal(b1-4)GlcNAc(b1-"
  for (encode in list(fp_sweetnet, fp_glycangt, fp_gifflar)) {
    for (opts in list(
      list(max_size = 2),
      list(unknown = "error"),
      list(method = "sweetnet")
    )) {
      expect_error(
        do.call(encode, c(list(x, checkpoint = "local.pt"), opts)),
        "must be empty"
      )
    }
    for (path in list(NULL, NA_character_, character(), "")) {
      expect_error(encode(x, checkpoint = path), class = "glyfp_error_argument")
    }
    expect_error(
      encode(NA_character_, checkpoint = "local.pt", na_action = "error"),
      class = "glyfp_error_missing"
    )
    expect_error(
      encode(x, checkpoint = "local.pt", na_action = "drop"),
      class = "glyfp_error_argument"
    )
  }
  expect_error(fp_sweetnet(x, download = NA), class = "glyfp_error_argument")
  expect_error(fp_sweetnet(x, cache_dir = NULL), class = "glyfp_error_argument")
  for (encode in list(fp_glycangt)) {
    expect_error(encode(x), class = "glyfp_error_argument")
    expect_error(
      encode(x, checkpoint = tempfile()),
      class = "glyfp_error_backend"
    )
    expect_error(
      encode(x, checkpoint = "https://example.org/model.pt"),
      class = "glyfp_error_argument"
    )
    for (opts in list(
      list(download = TRUE),
      list(download = FALSE),
      list(cache_dir = tempfile())
    )) {
      expect_error(
        do.call(encode, c(list(x, checkpoint = "local.pt"), opts)),
        "must be empty"
      )
    }
  }
  for (encode in list(fp_motif_binary, fp_motif_count)) {
    expect_error(encode(x, checkpoint = "model.pt"), "must be empty")
    expect_error(encode(x, download = FALSE), "must be empty")
    expect_error(encode(x, cache_dir = tempfile()), "must be empty")
  }
})

for (model in c("sweetnet-v1", "sweetnet-v2", "glycangt", "gifflar")) {
  local({
    model <- model
    test_that(
      paste(
        model,
        "matches official reference and independent batch contracts"
      ),
      {
        variable <- switch(
          model,
          `sweetnet-v1` = "GLYFP_TEST_CHECKPOINT",
          `sweetnet-v2` = "GLYFP_TEST_V2_CHECKPOINT",
          glycangt = "GLYFP_TEST_GLYCANGT_CHECKPOINT",
          gifflar = "GLYFP_TEST_GIFFLAR_CHECKPOINT"
        )
        checkpoint <- Sys.getenv(variable)
        skip_if(
          !nzchar(checkpoint),
          paste("Set", variable, "for real-weight validation")
        )
        skip_if_not_installed("torch")
        skip_if_not(
          torch::torch_is_installed(),
          "R torch runtime is not installed"
        )
        if (model == "gifflar") {
          reference <- readRDS(test_path("fixtures", "gifflar-reference.rds"))
          inputs <- vapply(reference, `[[`, "", "canonical")
          expected <- do.call(rbind, lapply(reference, `[[`, "embedding"))
        } else {
          file <- switch(
            model,
            `sweetnet-v1` = "sweetnet-reference.csv",
            `sweetnet-v2` = "sweetnet-v2-reference.csv",
            glycangt = "glycangt-reference.csv"
          )
          reference <- read.csv(
            test_path("fixtures", file),
            check.names = FALSE
          )
          inputs <- if (model == "glycangt") reference$key else reference$glycan
          expected <- as.matrix(reference[,
            if (model == "glycangt") -(1:2) else -1L
          ])
        }
        check_embedding_contract(
          if (startsWith(model, "sweetnet")) "sweetnet" else model,
          checkpoint,
          inputs,
          expected
        )
      }
    )
  })
}
