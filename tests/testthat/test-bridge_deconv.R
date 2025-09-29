test_that("BRIDGEdeconv works with PAM50 reference", {

  data(reference_BRIDGE_PAM50, package = utils::packageName())
  # create synthetic expr matrix with sufficient overlap
  set.seed(3)
  gsel <- sample(rownames(reference_BRIDGE_PAM50), 600)
  expr <- matrix(round(runif(600 * 3), 3)*1000, nrow = 600,
                 dimnames = list(gsel, paste0("P", 1:3)))

  res <- BRIDGEdeconv(expr_matrix = expr, reference = "PAM50", bcor = FALSE)

  expect_true(is.matrix(res))
  expect_equal(nrow(res), ncol(expr))   # one row per sample
  expect_equal(colnames(res), colnames(reference_BRIDGE_PAM50)) # subtypes
  expect_true(all(is.finite(res)))
})

test_that("BRIDGEdeconv fails if expr_matrix is not a matrix", {
  expect_error(
    BRIDGEdeconv(expr_matrix = as.data.frame(matrix(1:10, nrow = 5)),
                 reference = "PAM50"),
    regexp = "must be a matrix"
  )
})

test_that("BRIDGEdeconv fails if bcor is not logical", {
  data(reference_BRIDGE_PAM50, package = utils::packageName())
  gsel <- sample(rownames(reference_BRIDGE_PAM50), 50)
  expr <- matrix(round(runif(50 * 2), 3)*1000, nrow = 50,
                 dimnames = list(gsel, c("S1","S2")))

  expect_error(
    BRIDGEdeconv(expr_matrix = expr, reference = "PAM50", bcor = "yes"),
    regexp = "`bcor` must be a single logical"
  )
})

test_that("BRIDGEdeconv errors if no overlap", {
  data(reference_BRIDGE_PAM50, package = utils::packageName())

  fake_genes <- paste0("FAKE", seq_len(20))
  expr <- matrix(round(runif(20 * 2), 3)*1000, nrow = 20,
                 dimnames = list(fake_genes, c("S1","S2")))

  expect_error(
    BRIDGEdeconv(expr_matrix = expr, reference = "PAM50", bcor = FALSE),
    regexp = "No overlapping genes"
  )
})

test_that("BRIDGEdeconv warns if overlap < 50%", {
  data(reference_BRIDGE_PAM50, package = utils::packageName())

  # take 30 genes from reference, add 70 fakes → 30/100 = 30% overlap
  overlap_genes <- head(rownames(reference_BRIDGE_PAM50), 30)
  fake_genes <- paste0("FAKE", seq_len(70))
  gsel <- c(overlap_genes, fake_genes)

  expr <- matrix(rnorm(length(gsel) * 2), nrow = length(gsel),
                 dimnames = list(gsel, c("S1","S2")))

  expect_warning(
    BRIDGEdeconv(expr_matrix = expr, reference = "PAM50", bcor = FALSE),
    regexp = "Low overlap"
  )
})
