test_that("GWAS plotting constructors return reusable plot objects", {
  gwas <- data.frame(
    chr = rep(c("1", "2"), each = 6),
    pos = rep(seq(100, 600, by = 100), 2),
    id = paste0("m", 1:12),
    p = c(0.5, 0.2, 0.03, 1e-5, 0.1, 0.7, 0.8, 0.4, 0.02, 2e-6, 0.06, 0.3),
    PVE = seq(0.01, 0.12, length.out = 12)
  )

  m <- plot_manhattan(gwas, point_size_by = "PVE", label_top = 2, draw = FALSE)
  q <- plot_qq(gwas, draw = FALSE)

  expect_s3_class(m, "manhattan_plot")
  expect_s3_class(q, "qq_plot")
})

test_that("regional plot keeps GWAS points aligned to LD variants", {
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
  gwas <- data.frame(
    chr = "1",
    pos = c(100, 200, 300, 400),
    id = paste0("s", 1:4),
    p = c(0.2, 1e-4, 0.03, 0.4)
  )

  p <- plot_ld_region(
    ld,
    metric = "r2",
    gwas = gwas,
    lead = "s2",
    cutline = 3,
    show_maf = FALSE,
    show_snp_connectors = TRUE,
    draw = FALSE
  )

  expect_s3_class(p, "ld_plot")
})
