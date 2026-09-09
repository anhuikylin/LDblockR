test_that("read_gwas parses headerless EMMAX .ps files", {
  tf <- tempfile(fileext = ".ps")
  writeLines(
    c(
      "chr1.S_3279 0.084 1.2e-08",
      "chr1.S_9451 -0.031 2.4e-05",
      "chr2:15000 0.012 0.004"
    ),
    tf
  )

  g <- read_gwas(tf)

  expect_equal(g$chr, c("chr1", "chr1", "chr2"))
  expect_equal(g$pos, c(3279, 9451, 15000))
  expect_equal(g$id, c("chr1.S_3279", "chr1.S_9451", "chr2:15000"))
  expect_equal(attr(g, "gwas_format"), "emmax")
  expect_true(all(g$p > 0 & g$p <= 1))
})

test_that("read_gwas recognizes GEMMA-style columns", {
  tab <- data.frame(
    chr = c(1, 1),
    rs = c("s1", "s2"),
    ps = c(100, 200),
    beta = c(0.1, -0.2),
    se = c(0.02, 0.03),
    p_wald = c(1e-6, 0.01),
    check.names = FALSE
  )

  g <- read_gwas(tab)

  expect_equal(g$pos, c(100, 200))
  expect_equal(g$id, c("s1", "s2"))
  expect_equal(g$beta, c(0.1, -0.2))
  expect_equal(g$se, c(0.02, 0.03))
  expect_equal(attr(g, "gwas_format"), "gemma")
})

test_that("read_gwas recognizes PLINK2 hash-prefixed headers", {
  tf <- tempfile(fileext = ".glm")
  writeLines(
    c(
      "## test metadata line",
      "#CHROM POS ID REF ALT A1 TEST OBS_CT BETA SE P",
      "1 123 v1 A G G ADD 100 0.2 0.04 1e-5",
      "2 456 v2 C T T ADD 100 -0.1 0.05 0.03"
    ),
    tf
  )

  g <- read_gwas(tf)

  expect_equal(g$chr, c("1", "2"))
  expect_equal(g$pos, c(123, 456))
  expect_equal(g$id, c("v1", "v2"))
  expect_equal(attr(g, "gwas_format"), "plink")
})
