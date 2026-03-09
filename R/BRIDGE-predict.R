#' BRIDGE Response Prediction
#'
#' Run deconvolution + predict therapy response scores for each patient,
#' given a bulk RNA-seq matrix and a specified clinical subtype.
#'
#' @param expr_matrix A genes × samples expression matrix (genes in rows, samples in columns).
#' @param subtype A character string indicating the breast cancer subtype.
#'   Must be one of: "ERpos_HER2neg", "HER2pos", "TNBC".
#' @param therapy A character string indicating the therapy type for response prediction.
#'   Must be one of: "CHEMO" for chemotherapy, "IMMUNO" for immunotherapy, "ANTI_HER2" for anti-her2, "ENDO" for endocrine therapy.
#' @param bcor Logical; whether to apply ComBat batch correction (default = TRUE).
#'
#' @return A list containing:
#'   \describe{
#'     \item{fractions}{A matrix of BRIDGE-estimated subtype abundances per sample.}
#'     \item{response_score}{A numeric vector of predicted response scores per sample.}
#'   }
#' @export
#'
BRIDGEpredict <- function(expr_matrix, subtype, therapy, bcor=T) {

  # check input
  .checkExprMatrix(expr_matrix)   # shared validation
  .checkSubtype(subtype)
  .checkTherapy(therapy)

  subtype <- match.arg(subtype, c("ERpos_HER2neg","HER2pos","TNBC","ERpos"))
  therapy <- match.arg(therapy, c("CHEMO","IMMUNO","ANTI_HER2","ENDO"))

  # select reference and model based on subtype
  if (subtype == "ERpos_HER2neg" & therapy == "CHEMO") {

    reference <- reference_BRIDGE_PAM50
    model     <- model_ERpos_HER2neg

  } else if (subtype == "ERpos_HER2neg" & therapy == "IMMUNO") {

    reference <- reference_BRIDGE_PAM50
    model     <- exploratory_model_ERpos_HER2neg_ICB

  } else if (subtype == "ERpos" & therapy == "ENDO") {

    reference <- reference_BRIDGE_PAM50
    model     <- exploratory_model_ERpos_ENDO

  } else if (subtype == "HER2pos" & therapy == "ANTI_HER2") {

    reference <- reference_BRIDGE_PAM50
    model     <- model_HER2pos

  } else if (subtype == "TNBC" & therapy == "CHEMO") {

    reference <- reference_BRIDGE_TNBC
    model     <- model_TNBC

  } else if (subtype == "TNBC" & therapy == "IMMUNO") {

    reference <- reference_BRIDGE_TNBC
    model     <- exploratory_model_TNBC_ICB

  } else {
    stop(
      "The input subtype and therapy combination have not been tested. ",
      "Please choose among:\n",
      "  subtype='ERpos_HER2neg' + therapy='CHEMO'\n",
      "  subtype='ERpos_HER2neg' + therapy='IMMUNO' [exploratory]\n",
      "  subtype='ERpos'         + therapy='ENDO'   [exploratory]\n",
      "  subtype='HER2pos'       + therapy='ANTI_HER2'\n",
      "  subtype='TNBC'          + therapy='CHEMO'\n",
      "  subtype='TNBC'          + therapy='IMMUNO' [exploratory]\n"
    )
  }

  # check overlap after reference is chosen
  .checkOverlap(expr_matrix, reference)

  # deconvolution
  fractions <- runDeconv(expr_matrix, bcor = bcor, sig=reference)

  # prediction
  response_score <- predict(model$mod, newdata = as.data.frame(fractions),
                            type='response')

  main_class <- ifelse(response_score >= model$thr, 'high', 'low')

  # predicted sub-classes
  sub_classes <- if (!is.null(model$qui)) {
    # extend breaks to -Inf/Inf to capture all scores outside training range
    breaks <- c(-Inf, model$qui[2:5], Inf)
    cut(response_score, breaks, labels = paste0('q', 1:5), include.lowest = TRUE)
  } else {
    rep(NA, length(response_score))
  }

  # return both
  list(
    fractions=fractions,
    BRIDGE_SCORE = data.frame(
      SCORE = response_score,
      CLASS = main_class,
      SUBCLASS = sub_classes
    )
  )
}
