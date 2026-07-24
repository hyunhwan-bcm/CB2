context("Regression extensions")

set.seed(20260723)

sample_data <- data.frame(
  dose = rep(seq(-1.5, 1.5, length.out = 6), 2),
  batch = factor(rep(c("a", "b"), each = 6))
)
total <- rep(50000, nrow(sample_data))
eta <- -7 + 0.8 * sample_data$dose +
  ifelse(sample_data$batch == "b", 0.15, 0)
rho <- 0.001
precision <- 1 / rho - 1
probability <- rbeta(
  nrow(sample_data),
  plogis(eta) * precision,
  (1 - plogis(eta)) * precision
)
count <- rbinom(nrow(sample_data), total, probability)

test_that("bbreg fits a continuous sample phenotype", {
  fit <- bbreg(count, total, ~ dose + batch, sample_data)
  expect_s3_class(fit, "bbreg")
  expect_true(fit$converged)
  expect_equal(fit$df.residual, nrow(sample_data) - 3)
  expect_true(all(c(
    "estimate", "std_error", "t_value", "df", "p_value"
  ) %in% colnames(fit$coefficient_table)))
  expect_true(is.finite(fit$coefficient_table["dose", "t_value"]))

  contrast <- bb_contrast(fit, c(dose = 1))
  expect_equal(contrast$estimate, unname(coef(fit)["dose"]))
  expect_equal(contrast$df, fit$df.residual)
})

test_that("Rcpp weighted products agree with base R", {
  x <- model.matrix(~ dose + batch, sample_data)
  weight <- runif(nrow(x), 0.1, 2)
  response <- rnorm(nrow(x))
  cpp <- getFromNamespace("bb_wls_system_cpp", "CB2")(
    x, weight, response
  )
  expect_equal(
    unname(cpp$information),
    unname(crossprod(x, weight * x)),
    tolerance = 1e-12
  )
  expect_equal(
    unname(drop(cpp$score_target)),
    unname(drop(crossprod(x, weight * response))),
    tolerance = 1e-12
  )
})

test_that("bb_screen preserves guide and gene annotation", {
  counts <- rbind(
    sg1 = count,
    sg2 = rbinom(length(total), total, 0.001),
    sg3 = rbinom(length(total), total, 0.002)
  )
  result <- bb_screen(
    counts, sample_data, ~ dose + batch, "dose",
    totals = total,
    gene = c("gene1", "gene1", "gene2")
  )
  expect_equal(nrow(result), 3)
  expect_equal(result$gene, c("gene1", "gene1", "gene2"))
  expect_true(all(result$fdr >= 0 & result$fdr <= 1, na.rm = TRUE))
})

test_that("negative-control calibration preserves raw inference", {
  control_t <- 1.8 * qt(
    (seq_len(100) - 0.5) / 100,
    df = 8
  )
  result <- data.frame(
    estimate = rep(1, 101),
    std_error = rep(1, 101),
    t_value = c(control_t, 5),
    df = rep(8, 101),
    p_value = 2 * pt(-abs(c(control_t, 5)), df = 8)
  )
  result$fdr <- p.adjust(result$p_value, method = "BH")
  calibrated <- bb_calibrate_controls(
    result,
    control = c(rep(TRUE, 100), FALSE)
  )
  expect_gt(attr(calibrated, "control_scale"), 1)
  expect_equal(calibrated$raw_t_value, result$t_value)
  expect_true(all(calibrated$p_value >= result$p_value))
  expect_lte(
    abs(mean(calibrated$p_value[seq_len(100)] < 0.05) - 0.05),
    0.02
  )
})

test_that("regression input errors are informative", {
  expect_error(
    bbreg(count, total, dose ~ batch, sample_data),
    "one-sided"
  )
  expect_error(
    bb_screen(matrix(count, nrow = 1), sample_data, ~ dose, "missing"),
    "one model-matrix coefficient"
  )
})
