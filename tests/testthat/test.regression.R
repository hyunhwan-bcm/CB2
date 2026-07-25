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

test_that("saturated GLS reproduces completed legacy CB2 summaries", {
  count_a <- matrix(c(74, 112, 91, 139), nrow = 1)
  total_a <- matrix(c(52000, 81000, 69000, 97000), nrow = 1)
  count_b <- matrix(c(128, 177, 164, 211), nrow = 1)
  total_b <- matrix(c(55000, 76000, 72000, 93000), nrow = 1)
  fit_a <- fit_ab(count_a, total_a)
  fit_b <- fit_ab(count_b, total_b)
  phat <- c(fit_a$phat, fit_b$phat)
  vhat <- c(fit_a$vhat, fit_b$vhat)

  x <- cbind(1, c(0, 1))
  gls <- getFromNamespace("bb_wls_solve_cpp", "CB2")(
    x, 1 / vhat, phat, covariance = TRUE
  )
  beta <- drop(gls$coefficient)
  covariance <- gls$covariance
  original_t <- diff(phat) / sqrt(sum(vhat))

  expect_equal(beta, c(phat[1], diff(phat)), tolerance = 1e-12)
  expect_equal(covariance[2, 2], sum(vhat), tolerance = 1e-12)
  expect_equal(
    beta[2] / sqrt(covariance[2, 2]),
    original_t,
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
  expect_true(all(c(
    "pearson_ratio", "scale", "dispersion_boundary", "converged"
  ) %in% names(result)))
  expect_type(result$dispersion_boundary, "logical")
  expect_true(all(is.finite(result$pearson_ratio[result$converged])))
  expect_true(all(is.finite(result$scale[result$converged])))
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

test_that("guide-consistency scores expose their exploratory empirical null", {
  consistency_input <- data.frame(
    gene = rep(c("control_a", "control_b", "null", "signal"), each = 5),
    estimate = c(
      -0.08, -0.04, 0, 0.04, 0.08,
      -0.10, -0.05, 0, 0.05, 0.10,
      -0.05, -0.02, 0.01, 0.04, 0.02,
      0.58, 0.62, 0.65, 0.68, 0.72
    ),
    std_error = rep(0.10, 20),
    converged = TRUE
  )
  control <- consistency_input$gene %in% c("control_a", "control_b")
  result <- bb_gene_consistency(
    consistency_input,
    control = control,
    min_control_genes = 2
  )

  expect_equal(nrow(result), 4)
  expect_true(all(c(
    "raw_statistic", "statistic", "guide_direction_agreement",
    "converged_fraction", "control_gene", "p_value", "fdr"
  ) %in% names(result)))
  expect_gte(attr(result, "null_scale"), 1)
  expect_equal(attr(result, "control_genes"), 2)
  expect_lt(
    result$p_value[result$gene == "signal"],
    result$p_value[result$gene == "null"]
  )
  expect_equal(
    result$guide_direction_agreement[result$gene == "signal"],
    1
  )
  expect_match(
    attr(result, "null_assumption"),
    "not biological-replicate inference",
    fixed = TRUE
  )
})

test_that("guide-consistency excludes failed fits and uses raw uncertainty", {
  result <- data.frame(
    gene = rep(c("control_a", "control_b", "signal"), each = 5),
    estimate = c(
      -0.08, -0.04, 0, 0.04, 0.08,
      -0.10, -0.05, 0, 0.05, 0.10,
      -10, 0.60, 0.60, 0.60, 0.60
    ),
    std_error = rep(0.30, 15),
    raw_std_error = rep(0.10, 15),
    converged = c(rep(TRUE, 10), FALSE, rep(TRUE, 4))
  )
  control <- result$gene %in% c("control_a", "control_b")
  gene_result <- bb_gene_consistency(
    result,
    control = control,
    min_control_genes = 2
  )

  signal <- gene_result[gene_result$gene == "signal", ]
  expect_equal(signal$n_guides, 4)
  expect_equal(signal$converged_fraction, 0.8)
  expect_equal(signal$estimate, 0.60, tolerance = 1e-12)
  expect_equal(signal$std_error, 0.05, tolerance = 1e-12)
})

test_that("nonconverged fits cannot produce inferential results", {
  fit <- bbreg(
    count, total, ~ dose, sample_data,
    maxit = 1L,
    tolerance = 1e-300
  )
  expect_false(fit$converged)
  expect_true(all(is.na(fit$coefficient_table[
    , c("std_error", "t_value", "p_value")
  ])))
  expect_true(all(is.na(vcov(fit))))

  contrast <- bb_contrast(fit, c(dose = 1))
  expect_true(is.finite(contrast$estimate))
  expect_true(all(is.na(contrast[
    , c("std_error", "t_value", "df", "p_value")
  ])))

  screen <- bb_screen(
    rbind(guide_a = count),
    totals = total,
    data = sample_data,
    formula = ~ dose,
    term = "dose",
    maxit = 1L,
    tolerance = 1e-300
  )
  expect_false(screen$converged)
  expect_true(all(is.na(screen[
    , c("estimate", "std_error", "t_value", "p_value", "fdr")
  ])))
  expect_true(all(c(
    "pearson_ratio", "scale", "dispersion_boundary"
  ) %in% names(screen)))
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
