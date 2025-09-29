test_that("BRIDGE_P runs and returns fractions + response_score", {
  # Choose the clinical subtype you know has a model in sysdata.rda.
  # Adjust if your internal models differ.
  subtype <- "ERpos_HER2neg"

  # We still need a reference-driven gene set to build expr.
  # Use PAM50 genes to synthesize inputs (deconv uses them internally anyway).
  data(reference_BRIDGE_PAM50, package = utils::packageName())

  set.seed(3)
  gsel <- sample(rownames(reference_BRIDGE_PAM50), 600)
  expr <- matrix(round(runif(600 * 3), 3)*1000, nrow = 600,
                 dimnames = list(gsel, paste0("P", 1:3)))

  # Smoke test: BRIDGEpredict uses internal models (sysdata.rda)
  res <- BRIDGEpredict(expr_matrix = expr, bcor=F, subtype = subtype)

  # structure
  expect_type(res, "list")
  expect_true(all(c("fractions", "BRIDGE_SCORE") %in% names(res)))

  frac  <- res$fractions
  score <- res$BRIDGE_SCORE

  testthat::expect_type(frac, "double")   # atomic type
  expect_true(is.matrix(frac))            # check matrix structure
  expect_equal(nrow(frac), ncol(expr))
  expect_true(is.numeric(score[,1]))
  expect_true(is.character(score[,2]))
  expect_equal(nrow(score), ncol(expr))

  # numeric sanity
  expect_true(all(is.finite(frac)))
  expect_true(all(is.finite(score[,1])))

})

