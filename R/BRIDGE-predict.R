#' BRIDGE Response Prediction
#'
#' Run deconvolution + predict therapy response scores for each patient,
#' given a bulk RNA-seq matrix and a specified clinical subtype.
#'
#' @param expr_matrix A genes × samples expression matrix (genes in rows, samples in columns).
#' @param subtype A character string indicating the breast cancer subtype.
#'   Must be one of: "ERpos_HER2neg", "HER2pos", "TNBC".
#' @param bcor Logical; whether to apply ComBat batch correction (default = TRUE).
#'
#' @return A list containing:
#'   \describe{
#'     \item{fractions}{A matrix of BRIDGE-estimated subtype abundances per sample.}
#'     \item{response_score}{A numeric vector of predicted response scores per sample.}
#'   }
#' @export
#'
BRIDGEpredict <- function(expr_matrix, subtype, bcor=T) {

  # check input
  .checkExprMatrix(expr_matrix)   # shared validation
  .checkSubtype(subtype)

  subtype <- match.arg(subtype, c("ERpos_HER2neg","HER2pos","TNBC"))

  # select reference and model based on subtype
  if (subtype == "ERpos_HER2neg") {
    reference <- reference_BRIDGE_PAM50
    model     <- model_ERpos_HER2neg
  } else if (subtype == "HER2pos") {
    reference <- reference_BRIDGE_PAM50
    model     <- model_HER2pos
  } else if (subtype == "TNBC") {
    reference <- reference_BRIDGE_TNBC
    model     <- model_TNBC
  }

  # check overlap after reference is chosen
  .checkOverlap(expr_matrix, reference)

  # deconvolution
  fractions <- runDeconv(expr_matrix, bcor = bcor, sig=reference)

  # prediction
  response_score <- predict(model$mod, newdata = as.data.frame(fractions),
                            type='response')

  # return both
  list(
    fractions=fractions,
    BRIDGE_SCORE = data.frame(
      SCORE = response_score,
      CLASS = ifelse(response_score >= model$thr, 'high', 'low')
    )
  )
}
