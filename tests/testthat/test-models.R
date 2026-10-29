test_that("model options and local-only checkpoints are enforced", {
  local_mocked_bindings(
    download.file = function(...) stop("Unexpected network"),
    .package = "utils"
  )
  x <- "Gal(b1-4)GlcNAc(b1-"
  for (method in c("sweetnet", "glycangt", "gifflar")) {
    for (opts in list(
      list(max_size = 2),
      list(unknown = "ignore"),
      list(unknown = "warn"),
      list(download = NA),
      list(cache_dir = NULL)
    )) {
      expect_error(
        do.call(
          glycan_fingerprint,
          c(list(x, method, checkpoint = "local.pt"), opts)
        ),
        class = "glyfp_error_argument"
      )
    }
    for (path in list(NULL, NA_character_, character(), "")) {
      expect_error(
        glycan_fingerprint(x, method, checkpoint = path),
        class = "glyfp_error_argument"
      )
    }
    expect_error(
      glycan_fingerprint(
        NA_character_,
        method,
        checkpoint = "local.pt",
        na_action = "error"
      ),
      class = "glyfp_error_missing"
    )
  }
  for (method in c("glycangt", "gifflar")) {
    expect_error(glycan_fingerprint(x, method), class = "glyfp_error_argument")
    expect_error(
      glycan_fingerprint(x, method, checkpoint = tempfile()),
      class = "glyfp_error_backend"
    )
    expect_error(
      glycan_fingerprint(
        x,
        method,
        checkpoint = "https://example.org/model.pt"
      ),
      class = "glyfp_error_argument"
    )
    expect_error(
      glycan_fingerprint(x, method, checkpoint = "local.pt", download = TRUE),
      class = "glyfp_error_argument"
    )
  }
  expect_error(
    glycan_fingerprint(x, checkpoint = "model.pt"),
    class = "glyfp_error_argument"
  )
  expect_false("python" %in% names(formals(glycan_fingerprint)))
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
