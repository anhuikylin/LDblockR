test_that("core LD statistics remain symmetric and bounded", {
  g <- matrix(
    c(
      0, 0, 1, 1, 2, 2,
      0, 1, 1, 1, 2, 2,
      2, 2, 1, 1, 0, 0,
      0, 0, 0, 1, 1, 2
    ),
    nrow = 6,
    ncol = 4
  )
  map <- data.frame(
    chr = rep("1", 4),
    pos = c(100, 200, 300, 400),
    id = paste0("s", 1:4)
  )
  x <- as_ld_data(g, map, sample_ids = paste0("sample", 1:6))
  ld <- ld_compute(x, measure = "both", min_n = 3)

  expect_s3_class(ld, "ld_result")
  expect_equal(ld$r2, t(ld$r2), tolerance = 1e-12)
  expect_equal(ld$dprime, t(ld$dprime), tolerance = 1e-12)
  expect_true(all(ld$r2[is.finite(ld$r2)] >= 0 & ld$r2[is.finite(ld$r2)] <= 1))
  expect_true(all(ld$dprime[is.finite(ld$dprime)] >= 0 & ld$dprime[is.finite(ld$dprime)] <= 1))
  expect_equal(diag(ld$r2), rep(1, 4))
  expect_equal(diag(ld$dprime), rep(1, 4))
})

test_that("fixed blocks and tag SNPs use the same variant coordinates", {
  g <- matrix(
    c(
      0, 0, 1, 1, 2, 2,
      0, 1, 1, 1, 2, 2,
      2, 2, 1, 1, 0, 0,
      0, 0, 0, 1, 1, 2
    ),
    nrow = 6,
    ncol = 4
  )
  map <- data.frame(chr = "1", pos = c(100, 200, 300, 400), id = paste0("s", 1:4))
  x <- as_ld_data(g, map, sample_ids = paste0("sample", 1:6))
  ld <- ld_compute(x, measure = "r2", min_n = 3)

  blocks <- detect_ld_blocks(
    ld,
    method = "fixed",
    fixed = data.frame(chr = "1", start = 100, end = 400),
    min_snps = 2
  )
  tags <- select_tag_snps(ld, threshold = 0.5, blocks = blocks)

  expect_s3_class(blocks, "ld_blocks")
  expect_equal(nrow(blocks), 1L)
  expect_true(nrow(tags) >= 1L)
  expect_true(all(tags$tag %in% x$variants$id))
})
