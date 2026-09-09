.parse_gff_attributes <- function(x) {
  parts <- strsplit(x, ";", fixed = TRUE)[[1L]]
  out <- character()
  for (part in parts) {
    part <- trimws(part)
    if (!nzchar(part)) next
    if (grepl("=", part, fixed = TRUE)) {
      z <- strsplit(part, "=", fixed = TRUE)[[1L]]
      key <- z[1L]
      value <- paste(z[-1L], collapse = "=")
    } else {
      z <- strsplit(part, "[[:space:]]+", perl = TRUE)[[1L]]
      key <- z[1L]
      value <- paste(z[-1L], collapse = " ")
    }
    value <- gsub('^"|"$', "", value)
    out[key] <- utils::URLdecode(value)
  }
  out
}

#' Read regional GFF3/GTF annotation
#'
#' @param file GFF3, GTF, or compressed annotation file.
#' @param region Optional genomic region.
#' @param feature_types Optional feature types to retain.
#' @return A normalized annotation data frame.
#' @export
read_gff3 <- function(file, region = NULL,
                      feature_types = c("gene", "mRNA", "transcript", "exon", "CDS",
                                        "five_prime_UTR", "three_prime_UTR", "UTR")) {
  reg <- parse_region(region)
  con <- .open_text(file)
  on.exit(close(con), add = TRUE)
  rows <- list()
  k <- 0L
  done <- FALSE
  while (!done) {
    lines <- readLines(con, n = 10000L, warn = FALSE)
    if (!length(lines)) break
    lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
    fields <- strsplit(lines, "\t", fixed = TRUE)
    for (f in fields) {
      if (length(f) < 9L) next
      chr <- f[1L]
      start <- suppressWarnings(as.numeric(f[4L]))
      end <- suppressWarnings(as.numeric(f[5L]))
      if (!is.finite(start) || !is.finite(end)) next
      if (!is.null(reg)) {
        if (chr != reg$chr[1L] || end < reg$start[1L]) next
        if (chr == reg$chr[1L] && start > reg$end[1L]) {
          done <- TRUE
          break
        }
      }
      if (!is.null(feature_types) && !f[3L] %in% feature_types) next
      a <- .parse_gff_attributes(f[9L])
      pick <- function(keys) {
        z <- a[keys]
        z <- z[!is.na(z) & nzchar(z)]
        if (length(z)) unname(z[1L]) else NA_character_
      }
      k <- k + 1L
      rows[[k]] <- data.frame(
        chr = chr, source = f[2L], type = f[3L], start = start, end = end,
        score = suppressWarnings(as.numeric(f[6L])), strand = f[7L], phase = f[8L],
        id = pick(c("ID", "transcript_id", "gene_id")),
        parent = pick(c("Parent", "gene_id")),
        name = pick(c("Name", "gene_name", "gene", "locus_tag", "ID")),
        attributes = f[9L], stringsAsFactors = FALSE
      )
    }
  }
  if (!k) {
    return(data.frame(chr = character(), source = character(), type = character(),
                      start = numeric(), end = numeric(), score = numeric(),
                      strand = character(), phase = character(), id = character(),
                      parent = character(), name = character(), attributes = character(),
                      stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, rows)
  out <- out[order(.chromosome_rank(out$chr), out$chr, out$start, out$end), , drop = FALSE]
  rownames(out) <- NULL
  class(out) <- c("ld_genes", class(out))
  out
}

.gwas_file_layout <- function(file) {
  con <- .open_text(file)
  on.exit(close(con), add = TRUE)
  skip <- 0L
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) return(list(header = FALSE, skip = skip, comment = "#"))
    txt <- trimws(line)
    if (!nzchar(txt)) {
      skip <- skip + 1L
      next
    }
    if (grepl("^#CHROM(?:[[:space:]]|$)", txt, ignore.case = TRUE)) {
      return(list(header = TRUE, skip = skip, comment = ""))
    }
    if (startsWith(txt, "#")) {
      skip <- skip + 1L
      next
    }
    z <- strsplit(txt, "[[:space:],]+", perl = TRUE)[[1L]]
    header <- any(grepl(
      "^(chr|chrom|chromosome|pos|position|bp|ps|p|pval|pvalue|p_value|p\\.value|p_wald|p_lrt|p_score|logp|neglog10p|snp|marker|rs|id)$",
      tolower(z)
    ))
    return(list(header = header, skip = skip, comment = "#"))
  }
}

.has_header <- function(file) .gwas_file_layout(file)$header

.gwas_guess_format <- function(tab, source = NULL) {
  low <- if (is.null(source)) "" else tolower(basename(source))
  nms <- tolower(names(tab))
  if (grepl("\\.ps(?:\\.(?:txt|gz))?$", low)) return("emmax")
  if (any(c("p_wald", "p_lrt", "p_score") %in% nms) || all(c("rs", "ps") %in% nms)) return("gemma")
  if (any(c("#chrom", "obs_ct", "test") %in% nms)) return("plink")
  if (all(c("snp", "chromosome", "position") %in% nms) && any(c("p.value", "p_value", "p") %in% nms)) return("gapit")
  if (all(c("marker", "chr", "pos") %in% nms) && any(c("trait", "f") %in% nms)) return("tassel")
  "generic"
}

.gwas_coords_from_id <- function(id) {
  z <- trimws(as.character(id))
  pos <- rep(NA_real_, length(z))
  chr <- rep(NA_character_, length(z))
  # Common forms include chr1.S_3279, chr1:3279, 1_3279 and SNP_chr1_3279.
  has_pos <- grepl("[._:][0-9]+(?:[._:][A-Za-z]+)*$", z, perl = TRUE)
  if (any(has_pos)) {
    pos_text <- sub("^.*[._:]([0-9]+)(?:[._:][A-Za-z]+)*$", "\\1", z[has_pos], perl = TRUE)
    pos[has_pos] <- suppressWarnings(as.numeric(pos_text))
    chr_text <- sub("^.*?((?:chr)?(?:[0-9]+|X|Y|Z|W|M|MT))[._:].*$", "\\1",
                    z[has_pos], perl = TRUE, ignore.case = TRUE)
    valid_chr <- grepl("^(?:chr)?(?:[0-9]+|X|Y|Z|W|M|MT)$", chr_text,
                       perl = TRUE, ignore.case = TRUE)
    idx <- which(has_pos)
    chr[idx[valid_chr]] <- chr_text[valid_chr]
  }
  data.frame(chr = chr, pos = pos, stringsAsFactors = FALSE)
}

#' Read regional association statistics
#'
#' @param x File or data frame.
#' @param region Optional region.
#' @param chr_col,pos_col,p_col,id_col Explicit column names; auto-detected by default.
#' @param value_is_logp Whether the selected value is already -log10(P).
#' @param header Whether a file has a header; auto-detected by default.
#' @param sep File separator; empty means arbitrary whitespace.
#' @return A normalized data frame with `chr`, `pos`, `p`, `logp`, and `id`.
#' @details A bare `regional_gwas.tsv` filename resolves to the bundled
#' regional example after installation. Paths that include a directory are
#' treated as user-supplied paths. EMMAX `.ps` files and common GAPIT, GEMMA,
#' PLINK/PLINK2, and TASSEL column names are recognized automatically. When
#' chromosome and position columns are absent, coordinates are recovered from
#' common marker-ID forms such as `chr1.S_3279`, `chr1:3279`, and `1_3279`.
#' @export
read_gwas <- function(x, region = NULL, chr_col = NULL, pos_col = NULL,
                      p_col = NULL, id_col = NULL, value_is_logp = FALSE,
                      header = NULL, sep = "") {
  source_path <- NULL
  if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x)) {
      bundled <- system.file("extdata", basename(x), package = "LDblockR")
      if (identical(dirname(x), ".") && nzchar(bundled) && file.exists(bundled)) {
        x <- bundled
      } else {
        .stopf("GWAS file does not exist: %s. Use example_data(\"regional\") for the built-in regional example.", x)
      }
    }
    source_path <- normalizePath(x, mustWork = FALSE)
    layout <- .gwas_file_layout(x)
    if (is.null(header)) header <- layout$header
    con <- .open_text(x)
    on.exit(close(con), add = TRUE)
    tab <- utils::read.table(
      con,
      header = header,
      sep = sep,
      quote = "",
      comment.char = layout$comment,
      skip = layout$skip,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    if (!header) {
      if (ncol(tab) < 3L) .stopf("Headerless GWAS input needs at least three columns.")
      if (grepl("\\.ps(?:\\.(?:txt|gz))?$", tolower(basename(source_path)))) {
        names(tab)[1:3] <- c("id", "effect", "p")
      } else {
        names(tab)[1:3] <- c("chr", "pos", "p")
      }
    }
  } else {
    tab <- as.data.frame(x, stringsAsFactors = FALSE)
  }
  if (!nrow(tab)) .stopf("GWAS table is empty.")

  detected_format <- .gwas_guess_format(tab, source_path)
  if (is.null(id_col)) {
    id_col <- .match_column(tab, c("id", "snp", "marker", "rsid", "rs", "name", "variant", "variant_id"), FALSE)
  }
  if (is.null(chr_col)) {
    chr_col <- .match_column(tab, c("chr", "chrom", "chromosome", "#chrom"), FALSE)
  }
  if (is.null(pos_col)) {
    pos_col <- .match_column(tab, c("pos", "position", "bp", "site", "ps", "base_pair", "basepair"), FALSE)
  }

  if ((is.null(chr_col) || is.null(pos_col)) && !is.null(id_col)) {
    coords <- .gwas_coords_from_id(tab[[id_col]])
    if (is.null(chr_col) && any(!is.na(coords$chr))) {
      tab$.LDblockR_chr <- coords$chr
      chr_col <- ".LDblockR_chr"
    }
    if (is.null(pos_col) && any(is.finite(coords$pos))) {
      tab$.LDblockR_pos <- coords$pos
      pos_col <- ".LDblockR_pos"
    }
  }
  if (is.null(chr_col)) {
    .stopf("Cannot find chromosome information. Supply chr_col or use marker IDs such as chr1.S_3279 or chr1:3279.")
  }
  if (is.null(pos_col)) {
    .stopf("Cannot find position information. Supply pos_col or use marker IDs containing genomic coordinates.")
  }

  if (is.null(p_col)) {
    p_col <- .match_column(
      tab,
      c("p", "pval", "pvalue", "p_value", "p.value", "p_wald", "p_lrt", "p_score",
        "p-value", "p_value_wald", "logp", "neglog10p", "minus_log10_p", "value"),
      TRUE, "P/value"
    )
    if (tolower(p_col) %in% c("logp", "neglog10p", "minus_log10_p")) value_is_logp <- TRUE
  }
  pve_col <- .match_column(tab, c("PVE", "pve", "var_explained", "variance_explained",
                                 "percent_variance_explained"), FALSE)
  beta_col <- .match_column(tab, c("beta", "effect", "estimate", "coefficient", "effect_size"), FALSE)
  se_col <- .match_column(tab, c("se", "stderr", "standard_error", "std_err"), FALSE)

  value <- suppressWarnings(as.numeric(as.character(tab[[p_col]])))
  if (value_is_logp) {
    logp <- value
    p <- 10^(-value)
  } else {
    p <- value
    logp <- .safe_log10(p)
  }
  out <- data.frame(
    chr = as.character(tab[[chr_col]]),
    pos = suppressWarnings(as.numeric(as.character(tab[[pos_col]]))),
    p = p,
    logp = logp,
    id = if (is.null(id_col)) paste0(tab[[chr_col]], ":", tab[[pos_col]]) else as.character(tab[[id_col]]),
    stringsAsFactors = FALSE
  )
  if (!is.null(pve_col)) out$PVE <- suppressWarnings(as.numeric(as.character(tab[[pve_col]])))
  if (!is.null(beta_col)) out$beta <- suppressWarnings(as.numeric(as.character(tab[[beta_col]])))
  if (!is.null(se_col)) out$se <- suppressWarnings(as.numeric(as.character(tab[[se_col]])))
  valid <- is.finite(out$pos) & is.finite(out$logp) & !is.na(out$chr) & nzchar(out$chr)
  if (!value_is_logp) valid <- valid & is.finite(out$p) & out$p >= 0 & out$p <= 1
  out <- out[valid & .in_region(out$chr, out$pos, parse_region(region)), , drop = FALSE]
  if (!nrow(out)) .stopf("No finite GWAS rows remain after coordinate and P-value normalization.")
  out <- out[order(.chromosome_rank(out$chr), out$chr, out$pos, out$id), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "gwas_format") <- detected_format
  attr(out, "source") <- source_path
  out
}
