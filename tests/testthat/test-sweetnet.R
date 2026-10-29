test_that("SweetNet normalization only removes unknown inferred half-linkages", {
  x <- glyparse::auto_parse(c("Gal(b1-4)GlcNAc", "GlcNAc", "Neu5Ac"))
  inputs <- sweetnet_inputs(x)
  expect_identical(
    vapply(inputs, `[[`, "", "iupac"),
    c("Gal(b1-4)GlcNAc", "GlcNAc", "Neu5Ac")
  )
  expect_identical(
    sweetnet_inputs(glyparse::auto_parse(as.character(x))),
    inputs
  )
  expect_identical(inputs[[2]]$edges, list())
  expect_identical(as.character(inputs[[1]]$labels), c("Gal", "GlcNAc", "b1-4"))
  expect_identical(inputs[[1]]$edges, list(c(1, 2), c(2, 0)))
  for (text in c(
    "Gal(a1-",
    "Gal(b1-",
    "Gal(?2-",
    "Gal(??-",
    "Gal(?1-4)GlcNAc",
    "Gal(b1-?)GlcNAc"
  )) {
    expect_error(
      sweetnet_inputs(glyparse::auto_parse(text)),
      class = "glyfp_error_unsupported"
    )
  }
  expect_identical(
    sweetnet_inputs(glyparse::auto_parse("Gal(?1-")),
    sweetnet_inputs(glyparse::auto_parse("Gal"))
  )
})
