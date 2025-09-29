#' BRIDGE Deconvolution
#'
#' Estimate BRIDGE subtype fractions from bulk RNA-seq data,
#' given a bulk expression matrix and a specified reference set.
#'
#' @param expr_matrix A genes × samples expression matrix (genes in rows, samples in columns).
#' @param reference A character string indicating which BRIDGE reference to use.
#'   Must be one of: "PAM50", "TNBC", "INTCLUST".
#' @param bcor Logical; whether to apply batch correction (default = TRUE).
#'
#' @return A matrix of BRIDGE-estimated subtype abundances per sample.
#' @export
#'
BRIDGEdeconv <- function(expr_matrix, reference = c("PAM50", "TNBC", "INTCLUST"), bcor = TRUE) {

  .checkExprMatrix(expr_matrix)

  reference <- match.arg(reference)

  if (!is.logical(bcor) || length(bcor) != 1) {
    stop("`bcor` must be a single logical (TRUE/FALSE).")
  }

  # select reference matrix
  if (reference == "PAM50") {
    sig <- reference_BRIDGE_PAM50
  } else if (reference == "TNBC") {
    sig <- reference_BRIDGE_TNBC
  } else if (reference == "INTCLUST") {
    sig <- reference_BRIDGE_INTCLUST
  }

  # check overlap after reference is chosen
  .checkOverlap(expr_matrix, sig)

  # run deconvolution
  fractions <- runDeconv(expr_matrix, bcor = bcor, sig=sig)

  return(fractions)
}
