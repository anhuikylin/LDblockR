test_that("ld_data preserves provenance and analysis history", {
  g <- matrix(
    c(
      0, 0, 1, 1, 2, 2,
      0, 1, 1, 1, 2, 2,
      2, 2, 1, 1, 0, 0
    ),
    nrow = 6,
    ncol = 3
  )
  map <- data.frame(
    chr = c("1", "1", "1"),
    pos = c(100, 200, 300),
    id = c("s1", "s2", "s3")
  )

  x <- as_ld_data(g, map, sample_ids = paste0("sample", 1:6), source = "synthetic")

  expect_s3_class(x, "ld_data")
  expect_equal(x$provenance$source, "synthetic")
  expect_true(length(x$history) >= 1L)
  expect_equal(x$history[[1L]]$step, "construct")

  y <- filter_variants(
    x,
    min_maf = 0,
    max_maf = 0.5,
    max_missing = 1,
    max_het = 1,
    quiet = TRUE
  )
  expect_identical(y$provenance, x$provenance)
  expect_true(any(vapply(y$history, function(z) identical(z$step, "filter_variants"), logical(1L))))

  z <- subset_ld_data(y, variants = 1:2)
  expect_identical(z$provenance, x$provenance)
  expect_equal(z$n_variants, 2L)
  expect_true(any(vapply(z$history, function(q) identical(q$step, "subset_ld_data"), logical(1L))))
})
