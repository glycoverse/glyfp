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
