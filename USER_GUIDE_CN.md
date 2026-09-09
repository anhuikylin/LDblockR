# LDblockR 0.3.0 中文使用说明

LDblockR 是一个用于区域连锁不平衡（LD）分析、单倍型块识别、tagSNP 筛选、混合线性模型 GWAS 诊断和区域联合绘图的纯 R 软件包。

当前规范开发仓库：<https://github.com/anhuikylin/LDblockR>

要求 R >= 4.1.0。核心功能不依赖 Java、Bioconductor 或 ggplot2。

## 1. 安装

推荐使用 `pak`：

```r
pak::pak("anhuikylin/LDblockR")
library(LDblockR)
packageVersion("LDblockR")
```

## 2. 支持的数据

基因型输入：

- VCF / VCF.GZ
- HapMap
- PLINK BED/BIM/FAM
- PLINK PED/MAP
- 0/1/2/NA 数值矩阵

GWAS 结果可由 `read_gwas()` 统一读取，支持常见：

- EMMAX `.ps`
- GAPIT
- GEMMA
- PLINK / PLINK2
- TASSEL
- 通用 TSV/空格分隔表格

## 3. 内置示例

```r
regional <- example_data("regional")
regional

paths <- example_data("tassel")
paths
```

`regional` 包含区域 VCF、HapMap、GWAS、GFF3、固定 block、样本组、矩阵/map 和 PLINK 示例。

## 4. 读取区域基因型

### VCF

```r
x <- read_vcf_region(
  regional[["vcf"]],
  region = "chr1:1000000-1100000",
  min_maf = 0.01,
  max_missing = 0.25,
  quiet = TRUE
)
```

如果系统已安装 `bcftools` 且 VCF 已建立索引，区域提取会自动加速。

### HapMap

```r
x <- read_hapmap(
  regional[["hapmap"]],
  region = "chr1:1000000-1100000",
  min_maf = 0.01,
  quiet = TRUE
)
```

### PLINK

```r
x <- read_plink(
  regional[["plink"]],
  region = "chr1:1000000-1100000",
  min_maf = 0.01,
  quiet = TRUE
)
```

### 自动识别

```r
x <- read_genotypes(
  regional[["vcf"]],
  region = "chr1:1000000-1100000",
  min_maf = 0.01
)
```

## 5. `ld_data` 对象

0.3.0 开始，`ld_data` 会保存轻量的来源追踪和分析历史：

```r
x$provenance
x$history
```

过滤和子集操作会继续保留原始 provenance，并追加 history：

```r
x2 <- filter_variants(x, min_maf = 0.05)
x2$history

x3 <- subset_ld_data(x2, variants = 1:20)
x3$history
```

这样可以追踪一个区域对象是从什么文件、经过什么过滤和子集步骤得到的。

## 6. 计算 LD

```r
ld <- ld_compute(
  x,
  measure = "both",
  r2_method = "auto",
  min_n = 5
)
```

支持：

- dosage `r²`
- phased-haplotype `r²`
- EM-estimated `D′`
- 可选 `D′` profile-likelihood CI

如果只需要一定物理距离内的 LD：

```r
ld <- ld_compute(
  x,
  measure = "r2",
  max_distance = 200000
)
```

对于非常密集的区域，建议在读取阶段合理设置 `region`、`max_variants` 和 `on_excess`，避免生成不必要的大型 dense LD 矩阵。

## 7. 单倍型块

```r
blocks <- detect_ld_blocks(
  ld,
  method = "gabriel"
)
```

支持：

```text
gabriel
solid_spine
strong
four_gamete
fixed
none
```

固定坐标示例：

```r
fixed <- data.frame(
  chr = "1",
  start = 1000000,
  end = 1050000
)

blocks <- detect_ld_blocks(
  ld,
  method = "fixed",
  fixed = fixed
)
```

## 8. tagSNP

```r
tags <- select_tag_snps(
  ld,
  threshold = 0.8,
  blocks = blocks
)

tags
```

当提供 `blocks` 时，在各 block 内进行 greedy LD coverage 筛选。

## 9. 读取 GWAS 结果

### 常规表格

```r
gwas <- read_gwas("gwas.tsv")
```

只要存在常见的 chromosome、position、P-value 和 SNP ID 列，通常无需手动指定列名。

### EMMAX `.ps`

例如：

```text
chr1.S_3279  0.084  1.2e-08
chr1.S_9451 -0.031  2.4e-05
```

可直接：

```r
gwas <- read_gwas("shoot_tzr_salt.ps")
```

程序会从 `chr1.S_3279` 自动解析：

```text
chr = chr1
pos = 3279
```

也支持 `chr1:3279` 和 `1_3279` 等常见 ID。

### GEMMA

`chr / rs / ps / beta / se / p_wald` 等列可自动识别。

### PLINK2

`#CHROM / POS / ID / BETA / SE / P` 可直接读取，包括 `#CHROM` 开头的 header。

## 10. Manhattan 和 Q-Q 图

```r
m <- plot_manhattan(
  gwas,
  point_size_by = "PVE",
  label_top = 5,
  draw = FALSE
)

q <- plot_qq(
  gwas,
  draw = FALSE
)

save_gwas_plot(m, "Manhattan.pdf", width = 8.5, height = 5.4)
save_gwas_plot(q, "QQ.pdf", width = 6.5, height = 6.0)
```

Manhattan 图支持：

- 固定染色体间隔
- reference chromosome lengths
- chromosome/effect/PVE 着色
- PVE 点大小
- 多条 threshold/suggestive line
- top SNP 标注

## 11. 区域 GWAS + LD 联合图

```r
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
  square_cells = TRUE,
  show_snp_connectors = TRUE,
  heatmap_colors = c("#FFFDF2", "#FDBB55", "#B40426"),
  draw = FALSE
)

save_ld_plot(
  p,
  "regional_LD.pdf",
  width = 5.5,
  height = 7
)
```

区域图的 GWAS 点、SNP track、虚线 connector、block 和 LD heatmap 使用同一套 SNP index 对齐，避免上方 GWAS 与下方热图错位。

## 12. 混合线性模型 GWAS

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
  max_missing = 0.20,
  chunk_size = 256
)
```

模型采用 P3D-style single-variance-component MLM：先估计 null model 的方差比，再使用固定方差进行 SNP 扫描。

结果中：

- `PVE`：marker-level partial variance explained
- `model_PVE`：null model variance-component estimate

## 13. 从 GWAS 自动选区域

```r
region <- gwas_ld_region(
  gwas,
  cutline = 5.7,
  flank = 200000
)

region$region
region$lead_id
region$selected
```

然后用返回的区域读取 VCF/HapMap/PLINK 并绘制局部 LD。

## 14. 保存结果

```r
export_ld(ld, prefix = "my_region")
save_ld_plot(p, "regional_LD.svg", width = 5.5, height = 7)
save_ld_plot(p, "regional_LD.png", width = 5.5, height = 7, dpi = 600)
```

## 15. 验证与开发

仓库中的验证文件：

```text
inst/VALIDATION.md
inst/GWAS_VALIDATION.md
scripts/validate_ldblockshow.sh
scripts/reproduce_manuscript.R
```

核心 LD 结果已与官方 LDBlockShow 示例进行逐 pair 数值比较。

0.3.0 新增：

```text
tests/testthat/
.github/workflows/R-CMD-check.yaml
```

每次 push 到 `master` 后，会在 Linux、Windows 和 macOS 上自动执行 `R CMD check`。

## 16. 推荐分析顺序

```text
read_genotypes()
      ↓
ld_data
      ↓
filter_variants()
      ↓
ld_compute()
      ↓
detect_ld_blocks()
      ↓
select_tag_snps()
      ↓
read_gwas()
      ↓
plot_ld_region()
      ↓
save_ld_plot()
```

这条流程适合作为 LDblockR 的标准区域 GWAS/LD 分析流程。
