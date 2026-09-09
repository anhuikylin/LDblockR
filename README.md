# LDblockR

**LDblockR** is an R package for linkage-disequilibrium analysis, haplotype-block detection, tag-SNP selection, mixed-model GWAS diagnostics, and publication-ready regional association/LD visualization.

GitHub is the canonical development repository:
<https://github.com/anhuikylin/LDblockR>

## Highlights

- VCF/VCF.GZ, HapMap, PLINK BED/BIM/FAM, PLINK PED/MAP, and numeric genotype matrices.
- Dosage or phased-haplotype `r²` and EM-estimated `D′`.
- Gabriel, solid-spine, strong-pair, four-gamete, and fixed-coordinate haplotype blocks.
- Greedy tag-SNP selection, LD neighbors, LD decay, haplotype frequencies, subgroup comparisons, and batch regions.
- P3D-style single-variance-component mixed linear model GWAS with marker-level PVE and model-level PVE.
- Manhattan, Q-Q, and integrated regional GWAS + gene + block + LD heatmap figures.
- Vector output (PDF/SVG) and raster output (PNG/TIFF), with automatic rasterization for large heatmaps.
- Traceable `ld_data` objects containing source provenance and analysis history.
- Pure-R core with no Java or Bioconductor requirement.

## Installation

```r
pak::pak("anhuikylin/LDblockR")
library(LDblockR)
packageVersion("LDblockR")
```

The package requires R >= 4.1.0.

## 30-second regional example

```r
library(LDblockR)

regional <- example_data("regional")

x <- read_vcf_region(
  regional[["vcf"]],
  region = "chr1:1000000-1100000",
  min_maf = 0.01,
  quiet = TRUE
)

ld <- ld_compute(x, measure = "both")
blocks <- detect_ld_blocks(ld, method = "gabriel")
tags <- select_tag_snps(ld, threshold = 0.8, blocks = blocks)
gwas <- read_gwas(regional[["regional_gwas"]])

p <- plot_ld_region(
  ld,
  metric = "r2",
  gwas = gwas,
  genes = regional[["gff3"]],
  blocks = blocks,
  tags = tags,
  lead = gwas$id[which.max(gwas$logp)],
  cutline = 5,
  show_snp_connectors = TRUE,
  draw = FALSE
)

save_ld_plot(p, "regional_LD.pdf", width = 5.5, height = 7)
```

## GWAS result import

`read_gwas()` normalizes common association-result tables to `chr`, `pos`, `p`, `logp`, and `id`. It recognizes common GAPIT, GEMMA, PLINK/PLINK2, TASSEL, and EMMAX-style columns.

A headerless EMMAX `.ps` file such as:

```text
chr1.S_3279  0.084  1.2e-08
chr1.S_9451 -0.031  2.4e-05
```

can be read directly:

```r
gwas <- read_gwas("shoot_tzr_salt.ps")
head(gwas)
```

When chromosome and position columns are absent, marker IDs such as `chr1.S_3279`, `chr1:3279`, or `1_3279` are parsed automatically.

## Mixed-model GWAS

```r
paths <- example_data("tassel")
pheno <- read_tassel_phenotype(paths[["mdp_phenotype"]])

gwas <- gwas_mlm(
  genotype = paths[["mdp_genotype"]],
  phenotype = pheno,
  trait = "EarHT",
  covariates = c("location", "Q1", "Q2", "Q3"),
  replicate = "expand",
  min_maf = 0.05,
  max_missing = 0.20
)

m <- plot_manhattan(gwas, point_size_by = "PVE", label_top = 5, draw = FALSE)
q <- plot_qq(gwas, draw = FALSE)
```

## Traceable `ld_data`

Every imported/constructed `ld_data` object contains lightweight provenance and history records:

```r
x$provenance
x$history

x2 <- filter_variants(x, min_maf = 0.05)
x2$history
```

Filtering and subsetting preserve the original provenance and append a new history entry.

## Validation

The repository contains reproducible validation records and scripts under `inst/` and `scripts/`. The core phased example was numerically compared with official LDBlockShow output for all 861 variant pairs, with `r²` and `D′` agreeing to the precision reported by LDBlockShow.

The 0.3.0 development line adds automated `testthat` regression tests and GitHub Actions `R CMD check` across Linux, Windows, and macOS. The validated LD formulas are intentionally unchanged in this engineering release.

## Documentation

- Chinese user guide: [`USER_GUIDE_CN.md`](USER_GUIDE_CN.md)
- Core validation: [`inst/VALIDATION.md`](inst/VALIDATION.md)
- GWAS validation contract: [`inst/GWAS_VALIDATION.md`](inst/GWAS_VALIDATION.md)
- Reproducible scripts: [`scripts/`](scripts/)

## Development

Before committing changes to the numerical core or plotting coordinate system, add or update a regression test. In particular, preserve:

- SNP-to-heatmap one-to-one alignment;
- chromosome-width and fixed-gap behavior in Manhattan plots;
- block boundaries and tag-SNP coverage;
- symmetry/range constraints of `r²` and `D′`;
- reproducibility of the P3D-style GWAS scan.

## Citation

Run:

```r
citation("LDblockR")
```

for the package citation.

## License

MIT License.
