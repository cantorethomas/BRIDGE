utils::globalVariables(c(
  "reference_BRIDGE_PAM50",
  "reference_BRIDGE_TNBC",
  "reference_BRIDGE_INTCLUST"
))

#' Reference matrix for BRIDGE-PAM50
#'
#' Gene expression reference signatures used for BRIDGE deconvolution (PAM50).
#'
#' @format A numeric matrix with genes in rows and subtype features in columns.
#' @usage data(reference_BRIDGE_PAM50)
#' @examples
#' data(reference_BRIDGE_PAM50)
#' dim(reference_BRIDGE_PAM50)
"reference_BRIDGE_PAM50"

#' Reference matrix for BRIDGE-TNBC
#'
#' Gene expression reference signatures used for BRIDGE deconvolution (TNBC).
#'
#' @format A numeric matrix with genes in rows and subtype features in columns.
#' @usage data(reference_BRIDGE_TNBC)
#' @examples
#' data(reference_BRIDGE_TNBC)
#' dim(reference_BRIDGE_TNBC)
"reference_BRIDGE_TNBC"

#' Reference matrix for BRIDGE-INTCLUST
#'
#' Gene expression reference signatures used for BRIDGE deconvolution (IntClust).
#'
#' @format A numeric matrix with genes in rows and subtype features in columns.
#' @usage data(reference_BRIDGE_INTCLUST)
#' @examples
#' data(reference_BRIDGE_INTCLUST)
#' dim(reference_BRIDGE_INTCLUST)
"reference_BRIDGE_INTCLUST"

