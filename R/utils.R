#' Internal helper functions for BRIDGE
#'
#' These functions implement SVR-based deconvolution and utility checks.
#' They are not exported and are intended for internal package use only.
#'
#' @importFrom e1071 svm
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom stats sd predict
#' @noRd

.checkExprMatrix <- function(expr_matrix) {
  if (!is.matrix(expr_matrix)) {
    stop("`expr_matrix` must be a matrix with genes in rows and samples in columns.")
  }
  if (is.null(rownames(expr_matrix))) {
    stop("`expr_matrix` must have gene identifiers as row names.")
  }
  if (is.null(colnames(expr_matrix))) {
    stop("`expr_matrix` must have sample identifiers as column names.")
  }
  invisible(TRUE)
}

.checkSubtype <- function(subtype) {
  valid <- c("ERpos_HER2neg", "HER2pos", "TNBC")
  if (length(subtype) != 1 || !subtype %in% valid) {
    stop("`subtype` must be one of: ", paste(valid, collapse = ", "), ".")
  }
  invisible(TRUE)
}

.checkOverlap <- function(expr_matrix, reference, min_frac = 0.5) {
  overlap <- length(intersect(rownames(expr_matrix), rownames(reference)))
  total   <- nrow(reference)
  frac    <- overlap / total

  if (overlap == 0) {
    stop("No overlapping genes between expr_matrix and reference. Cannot run deconvolution.")
  }

  if (frac < min_frac) {
    warning(
      "Low overlap between expr_matrix and reference (",
      overlap, " out of ", total, ", ",
      round(frac * 100, 1), "%). ",
      "Results may be unreliable (recommended more than",
      min_frac * 100, "% overlap)."
    )
  }
  invisible(TRUE)
}


.runBCor <- function(dat, sig){

  cmg <- intersect(rownames(dat), rownames(sig))
  dat <- dat[cmg,]
  sig <- sig[cmg,]

  bc <- sva::ComBat(cbind(dat, sig), batch =
                      rep.int(c('a', 'b'),
                              c(ncol(dat), ncol(sig))))

  list(
    bc[,1:ncol(dat)],
    bc[,-c(1:ncol(dat))]
  )

}

.filterNA <- function(dat, sig){
  zeros <- apply(dat, 1, function(x) any(is.na(x)))
  kp <- which(!zeros)
  sig <- sig[kp,]
  dat <- dat[kp,]
  return(list(dat, sig))
}

.estimate_cell_frac_svm <- function(bulkp, mat_sig){
  # init check
  cell_types <- colnames(mat_sig)

  if (sd(bulkp) == 0){
    stop("The bulk mixture has no variations \n")
  }
  bulkp <- (bulkp - mean(bulkp)) / sd(bulkp)

  # model cell fraction using nu-SVM
  nu_list   <- c(0.25, 0.5, 0.75)
  frac_list <- array(0, dim = c(length(nu_list), length(cell_types)))
  err_list  <- rep(Inf, length(nu_list))
  for (i in 1:length(nu_list)){
    model_nu <- svm(bulkp ~ mat_sig, kernel = "linear",
                    type = "nu-regression", scale = F,
                    nu = nu_list[i])

    # adjust fractions
    frac_nu <- t(model_nu$coefs) %*% model_nu$SV
    frac_nu[which(frac_nu < 0)] <- 0
    if (sum(frac_nu) != 0){
      frac_nu <- frac_nu / sum(frac_nu)
    }

    # estimate model errors
    err_nu  <- sqrt(sum((bulkp - mat_sig %*% t(frac_nu))^2))

    frac_list[i, ] <- frac_nu
    err_list[i]    <- err_nu
  }

  # pick fractions with the least error
  frac <- frac_list[which.min(err_list), ]

  return(frac)
}


.runSVM <- function(bulk, signature, qn=F){

  # prepare data
  genes      <- intersect(rownames(signature), rownames(bulk))
  cell_types <- colnames(signature)
  samples    <- colnames(bulk)

  bulk <- bulk[genes, ]
  if (qn){
    bulk <- normalize.quantiles(bulk)
  }
  signature <- signature[genes, ]
  signature <- (signature - mean(signature)) / sd(as.vector(signature))
  # estimate cell fractions
  cell_frac <- array(0, dim = c(length(samples), length(cell_types)),
                     dimnames = list(samples, cell_types))
  for (i in 1:length(samples)){
    frac_i <- .estimate_cell_frac_svm(bulkp   = as.matrix(bulk[, i]),
                                      mat_sig = as.matrix(signature))
    cell_frac[i, ] <- as.matrix(frac_i)
  }

  return(cell_frac)

}

#' Internal: Run BRIDGE deconvolution engine
#'
#' Applies batch correction, NA filtering, and SVR-based deconvolution
#' using a given reference signature.
#'
#' @param dat A genes x samples expression matrix
#' @param bcor Logical; whether to apply ComBat batch correction
#' @param sig A reference signature matrix
#'
#' @return A matrix of estimated fractions (samples x subtypes)
#' @noRd
runDeconv <- function(dat, bcor=T, sig){

  sig <- log2(sig+1)
  if (any(apply(dat, 2, max) > 50)) dat <- log2(dat+1)

  if (bcor){
    # 3. Run Batch normalization
    all_dat <- .runBCor(dat, sig)
    sig <- all_dat[[2]]
    dat <- all_dat[[1]]
  } else{
    cmn <- intersect(rownames(dat), rownames(sig))
    sig <- sig[cmn,]
    dat <- dat[cmn,]
  }

  # 2. Filter NAs
  all_dat <- .filterNA(dat, sig)
  sig <- all_dat[[2]]
  dat <- all_dat[[1]]

  # 4. get Fractions
  .runSVM(dat, sig)

}

