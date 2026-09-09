# LDblockR 0.3.0

* Promotes GitHub as the canonical development repository and synchronizes package metadata with the current development line.
* Adds `testthat` regression tests for `ld_data`, core LD calculations, block/tag-SNP workflows, plotting objects, and common GWAS input styles.
* Adds GitHub Actions `R CMD check` on Linux, Windows, and macOS for every push and pull request.
* Extends `ld_data` with lightweight provenance and analysis-history records that survive filtering and subsetting.
* Improves `read_gwas()` auto-detection for EMMAX `.ps` output and common GAPIT, GEMMA, PLINK/PLINK2, and TASSEL column names. Coordinates can be recovered from marker IDs such as `chr1.S_3279`, `chr1:3279`, and `1_3279` when explicit chromosome/position columns are absent.
* Replaces the hard-coded startup version with the installed package version.
* Keeps the validated LD and P3D-style MLM numerical core unchanged while establishing regression infrastructure for future performance refactoring.

# LDblockR 0.0.1

* First public release for reproducible regional linkage-disequilibrium
  analysis in R.
* Supports VCF/VCF.GZ, HapMap, PLINK BED/BIM/FAM and PED/MAP, and numeric
  genotype matrices.
* Provides dosage or phased-haplotype r-squared, EM-estimated D-prime,
  configurable haplotype-block definitions, greedy tag-SNP selection, LD
  decay, haplotype frequencies, subgroup comparisons, and batch regional
  analysis.
* Provides integrated GWAS, GFF3/GTF, MAF, block, tag-SNP, and LD heatmap
  graphics with PDF, SVG, PNG, and TIFF export.
* Includes a one-call R workflow, a command-line interface, bundled example
  data, a reproduction script, and validation records.
