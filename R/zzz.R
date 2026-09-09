#' LDblockR package
#'
#' Regional linkage disequilibrium analysis and visualization.
#' @keywords internal
"_PACKAGE"

.onAttach <- function(libname, pkgname) {
  version <- tryCatch(as.character(utils::packageVersion(pkgname)), error = function(e) "development")
  packageStartupMessage(
    sprintf("LDblockR %s: regional LD, integrated visualization and GWAS diagnostics", version)
  )
}
