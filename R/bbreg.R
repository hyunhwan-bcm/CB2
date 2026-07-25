# Beta-binomial regression with t-based coefficient tests
#
# This is a dependency-free research implementation.  It fits the mean model
#
#   logit(mu_i) = x_i' beta
#
# using the beta-binomial variance
#
#   Var(K_i) = n_i mu_i (1 - mu_i) {1 + (n_i - 1) rho}.
#
# The coefficient estimates are obtained by feasible IRLS.  rho is estimated
# from the Pearson estimating equation, and coefficient tests use a Student
# t reference distribution with residual degrees of freedom.

#' @importFrom stats binomial glm.control glm.fit model.frame model.matrix
#'   na.fail plogis printCoefmat qlogis setNames terms uniroot
NULL

.bb_stop <- function(message) {
  stop(message, call. = FALSE)
}

.bb_validate_response <- function(count, total) {
  if (!is.numeric(count) || !is.numeric(total)) {
    .bb_stop("`count` and `total` must be numeric.")
  }
  if (length(count) != length(total) || length(count) < 2L) {
    .bb_stop("`count` and `total` must have the same length (at least two).")
  }
  if (anyNA(count) || anyNA(total) ||
      any(!is.finite(count)) || any(!is.finite(total))) {
    .bb_stop("`count` and `total` cannot contain missing or non-finite values.")
  }
  if (any(total <= 0) || any(count < 0) || any(count > total)) {
    .bb_stop("Each observation must satisfy 0 <= `count` <= `total`, with `total` > 0.")
  }
  integer_like <- function(x) all(abs(x - round(x)) < sqrt(.Machine$double.eps))
  if (!integer_like(count) || !integer_like(total)) {
    .bb_stop("`count` and `total` must contain integer-valued counts.")
  }
  invisible(TRUE)
}

.bb_make_design <- function(formula, data, n) {
  if (!inherits(formula, "formula") || length(formula) != 2L) {
    .bb_stop("`formula` must be a one-sided formula, for example `~ dose + batch`.")
  }
  if (!is.data.frame(data) || nrow(data) != n) {
    .bb_stop("`data` must be a data frame with one row per observation.")
  }
  mf <- model.frame(formula, data = data, na.action = na.fail)
  x <- model.matrix(formula, data = mf)
  qr_x <- qr(x)
  if (qr_x$rank < ncol(x)) {
    .bb_stop("The design matrix is not full rank; remove aliased covariates.")
  }
  if (nrow(x) <= ncol(x)) {
    .bb_stop("The model needs more observations than fitted coefficients.")
  }
  list(x = x, terms = terms(mf), contrasts = attr(x, "contrasts"))
}

.bb_inverse <- function(a) {
  chol2inv(chol(a))
}

.bb_wls_system <- function(x, weight, response) {
  if (exists("bb_wls_system_cpp", mode = "function", inherits = TRUE)) {
    return(bb_wls_system_cpp(x, weight, response))
  }
  list(
    information = crossprod(x, weight * x),
    score_target = crossprod(x, weight * response)
  )
}

.bb_wls_solve <- function(x, weight, response, covariance = FALSE) {
  if (exists("bb_wls_solve_cpp", mode = "function", inherits = TRUE)) {
    result <- bb_wls_solve_cpp(x, weight, response, covariance)
    result$coefficient <- drop(result$coefficient)
    return(result)
  }
  system <- .bb_wls_system(x, weight, response)
  coefficient <- drop(solve(system$information, system$score_target))
  if (!covariance) {
    return(list(coefficient = coefficient))
  }
  list(
    coefficient = coefficient,
    covariance = .bb_inverse(system$information)
  )
}

.bb_estimate_rho <- function(count, total, mu, df_residual,
                             upper = 1 - 1e-8) {
  binomial_variance <- total * mu * (1 - mu)
  pearson <- function(rho) {
    sum((count - total * mu)^2 /
          (binomial_variance * (1 + (total - 1) * rho)))
  }
  q0 <- pearson(0)
  if (!is.finite(q0)) {
    return(list(rho = 0, scale = 1, pearson = q0, boundary = TRUE))
  }
  if (q0 <= df_residual) {
    return(list(rho = 0, scale = 1, pearson = q0, boundary = FALSE))
  }
  qu <- pearson(upper)
  if (qu >= df_residual) {
    return(list(
      rho = upper,
      scale = max(1, qu / df_residual),
      pearson = qu,
      boundary = TRUE
    ))
  }
  rho <- uniroot(
    function(r) pearson(r) - df_residual,
    interval = c(0, upper),
    tol = 1e-10
  )$root
  list(rho = rho, scale = 1, pearson = pearson(rho), boundary = FALSE)
}

#' Fit beta-binomial regression for one guide
#'
#' @param count Guide counts, one value per sample.
#' @param total Total library counts, one value per sample.
#' @param formula One-sided mean-model formula, such as `~ dose + batch`.
#' @param data Sample-level covariates.
#' @param maxit Maximum feasible-IRLS iterations.
#' @param tolerance Relative convergence tolerance.
#' @param mu_bound Numerical bound applied to fitted proportions.
#' @return An object of class `bbreg`.
#' @export
bbreg <- function(count, total, formula, data, maxit = 100L,
                  tolerance = 1e-8, mu_bound = 1e-8) {
  .bb_validate_response(count, total)
  count <- as.numeric(count)
  total <- as.numeric(total)
  design <- .bb_make_design(formula, data, length(count))
  x <- design$x
  rank <- ncol(x)
  df_residual <- nrow(x) - rank

  initial <- suppressWarnings(glm.fit(
    x = x,
    y = count / total,
    weights = total,
    family = binomial(link = "logit"),
    control = glm.control(maxit = 50L, epsilon = tolerance)
  ))
  beta <- initial$coefficients
  if (any(!is.finite(beta))) {
    pooled <- (sum(count) + 0.5) / (sum(total) + 1)
    beta <- numeric(rank)
    beta[1L] <- qlogis(pooled)
  }

  converged <- FALSE
  rho_fit <- list(rho = 0, scale = 1, pearson = NA_real_,
                  boundary = FALSE)
  for (iteration in seq_len(maxit)) {
    eta <- drop(x %*% beta)
    mu <- pmin(pmax(plogis(eta), mu_bound), 1 - mu_bound)
    rho_fit <- .bb_estimate_rho(count, total, mu, df_residual)
    rho <- rho_fit$rho

    working_response <- eta + (count / total - mu) / (mu * (1 - mu))
    working_weight <- total * mu * (1 - mu) /
      (1 + (total - 1) * rho)
    beta_new <- tryCatch(
      .bb_wls_solve(
        x, working_weight, working_response, covariance = FALSE
      )$coefficient,
      error = function(e) rep(NA_real_, rank)
    )
    if (any(!is.finite(beta_new))) {
      .bb_stop("The IRLS update was singular; inspect sparse counts and the design.")
    }

    change <- max(abs(beta_new - beta) / pmax(1, abs(beta)))
    beta <- beta_new
    if (change < tolerance) {
      converged <- TRUE
      break
    }
  }

  eta <- drop(x %*% beta)
  mu <- pmin(pmax(plogis(eta), mu_bound), 1 - mu_bound)
  rho_fit <- .bb_estimate_rho(count, total, mu, df_residual)
  rho <- rho_fit$rho
  working_weight <- total * mu * (1 - mu) /
    (1 + (total - 1) * rho)
  final_wls <- .bb_wls_solve(
    x, working_weight, eta, covariance = TRUE
  )
  unscaled_covariance <- final_wls$covariance
  covariance <- rho_fit$scale * unscaled_covariance
  standard_error <- sqrt(diag(covariance))
  statistic <- beta / standard_error
  p_value <- 2 * pt(-abs(statistic), df = df_residual)
  if (!converged) {
    covariance[] <- NA_real_
    standard_error[] <- NA_real_
    statistic[] <- NA_real_
    p_value[] <- NA_real_
  }
  coefficient_table <- cbind(
    estimate = beta,
    std_error = standard_error,
    t_value = statistic,
    df = rep(df_residual, rank),
    p_value = p_value
  )
  rownames(coefficient_table) <- colnames(x)

  structure(list(
    coefficients = setNames(beta, colnames(x)),
    coefficient_table = coefficient_table,
    covariance = covariance,
    fitted.values = mu,
    linear.predictors = eta,
    residuals = count / total - mu,
    pearson = rho_fit$pearson,
    rho = rho,
    scale = rho_fit$scale,
    dispersion_boundary = rho_fit$boundary,
    df.residual = df_residual,
    rank = rank,
    count = count,
    total = total,
    formula = formula,
    data = data,
    design = x,
    weights = working_weight,
    converged = converged,
    iterations = iteration,
    call = match.call()
  ), class = "bbreg")
}

#' @export
coef.bbreg <- function(object, ...) {
  object$coefficients
}

#' @export
vcov.bbreg <- function(object, ...) {
  object$covariance
}

#' @export
fitted.bbreg <- function(object, ...) {
  object$fitted.values
}

#' @export
residuals.bbreg <- function(object, type = c("response", "pearson"), ...) {
  type <- match.arg(type)
  if (type == "response") {
    return(object$residuals)
  }
  variance <- object$fitted.values * (1 - object$fitted.values) *
    (1 + (object$total - 1) * object$rho) / object$total
  object$residuals / sqrt(variance)
}

#' @export
summary.bbreg <- function(object, ...) {
  ans <- list(
    call = object$call,
    coefficients = object$coefficient_table,
    rho = object$rho,
    pearson = object$pearson,
    df.residual = object$df.residual,
    converged = object$converged,
    iterations = object$iterations
  )
  class(ans) <- "summary.bbreg"
  ans
}

#' @export
print.summary.bbreg <- function(x, digits = max(3L, getOption("digits") - 3L),
                                ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  printCoefmat(x$coefficients, digits = digits, P.values = TRUE,
               has.Pvalue = TRUE)
  cat("\nBeta-binomial intraclass correlation (rho):",
      formatC(x$rho, digits = digits), "\n")
  cat("Pearson statistic / residual df:",
      formatC(x$pearson, digits = digits), "/",
      x$df.residual, "\n")
  cat("Converged:", x$converged, "after", x$iterations, "iterations\n")
  invisible(x)
}

#' Test a linear contrast of beta-binomial regression coefficients
#'
#' @param object A fitted `bbreg` object.
#' @param contrast Numeric contrast vector in coefficient order, or a named
#'   numeric vector whose names identify coefficients.
#' @param null Null value for the contrast.
#' @return A one-row data frame with estimate, standard error, t statistic,
#'   degrees of freedom, and two-sided p-value. If the underlying fit did not
#'   converge, the last-iterate contrast estimate is retained for diagnosis but
#'   its standard error, statistic, and p-value are missing.
#' @export
bb_contrast <- function(object, contrast, null = 0) {
  if (!inherits(object, "bbreg")) {
    .bb_stop("`object` must be a fitted `bbreg` object.")
  }
  p <- length(object$coefficients)
  if (!is.numeric(contrast) || anyNA(contrast)) {
    .bb_stop("`contrast` must be a numeric vector without missing values.")
  }
  if (!is.null(names(contrast))) {
    if (!all(names(contrast) %in% names(object$coefficients))) {
      .bb_stop("A named contrast contains an unknown coefficient.")
    }
    complete <- setNames(numeric(p), names(object$coefficients))
    complete[names(contrast)] <- contrast
    contrast <- complete
  }
  if (length(contrast) != p) {
    .bb_stop("An unnamed contrast must have one value per coefficient.")
  }
  estimate <- drop(crossprod(contrast, object$coefficients))
  if (!isTRUE(object$converged)) {
    return(data.frame(
      estimate = estimate,
      std_error = NA_real_,
      t_value = NA_real_,
      df = NA_real_,
      p_value = NA_real_,
      row.names = NULL
    ))
  }
  standard_error <- sqrt(drop(t(contrast) %*% object$covariance %*% contrast))
  statistic <- (estimate - null) / standard_error
  data.frame(
    estimate = estimate,
    std_error = standard_error,
    t_value = statistic,
    df = object$df.residual,
    p_value = 2 * pt(-abs(statistic), df = object$df.residual),
    row.names = NULL
  )
}

#' Apply beta-binomial regression guide by guide
#'
#' @param counts Guide-by-sample count matrix.
#' @param data Sample-level covariate data frame.
#' @param formula One-sided regression formula.
#' @param term Coefficient name to report.
#' @param totals Optional library-size vector. Defaults to column sums.
#' @param guide Optional guide identifiers. Defaults to row names.
#' @param gene Optional gene identifier per guide.
#' @param min_total_count Guides below this total count receive missing results.
#' @param ncores Number of forked workers on Unix-like systems. Windows uses
#'   one worker.
#' @param ... Additional arguments passed to `bbreg`.
#' @return A data frame with one row per guide, Benjamini--Hochberg FDR,
#'   convergence status, Pearson ratio, residual scale, and dispersion-boundary
#'   status. Nonconverged fits have missing inferential columns.
#' @export
bb_screen <- function(counts, data, formula, term, totals = NULL,
                      guide = rownames(counts), gene = NULL,
                      min_total_count = 10, ncores = 1L, ...) {
  if (!is.matrix(counts) && !is.data.frame(counts)) {
    .bb_stop("`counts` must be a numeric matrix or data frame.")
  }
  counts <- as.matrix(counts)
  storage.mode(counts) <- "double"
  if (anyNA(counts) || any(!is.finite(counts)) || any(counts < 0)) {
    .bb_stop("`counts` must contain finite, non-negative values.")
  }
  if (nrow(data) != ncol(counts)) {
    .bb_stop("`data` must have one row per count-matrix column.")
  }
  if (is.null(totals)) {
    totals <- colSums(counts)
  }
  if (length(totals) != ncol(counts)) {
    .bb_stop("`totals` must have one value per count-matrix column.")
  }
  if (any(counts > rep(totals, each = nrow(counts)))) {
    .bb_stop("A guide count cannot exceed its sample's `total`.")
  }
  if (is.null(guide)) {
    guide <- sprintf("guide_%d", seq_len(nrow(counts)))
  }
  if (length(guide) != nrow(counts) || anyDuplicated(guide)) {
    .bb_stop("`guide` must uniquely identify every row of `counts`.")
  }
  if (!is.null(gene) && length(gene) != nrow(counts)) {
    .bb_stop("`gene` must have one value per guide.")
  }
  if (length(ncores) != 1L || !is.finite(ncores) || ncores < 1) {
    .bb_stop("`ncores` must be one positive integer.")
  }
  ncores <- as.integer(ncores)

  design <- .bb_make_design(formula, data, ncol(counts))
  if (!term %in% colnames(design$x)) {
    .bb_stop(sprintf(
      "`term` must be one model-matrix coefficient: %s",
      paste(colnames(design$x), collapse = ", ")
    ))
  }

  one_guide <- function(i) {
    if (sum(counts[i, ]) < min_total_count) {
      return(c(estimate = NA_real_, std_error = NA_real_,
               t_value = NA_real_, df = NA_real_, p_value = NA_real_,
               rho = NA_real_, pearson_ratio = NA_real_, scale = NA_real_,
               dispersion_boundary = NA_real_,
               mean_cpm = mean(counts[i, ] / totals * 1e6),
               converged = 0))
    }
    fit <- tryCatch(
      bbreg(counts[i, ], totals, formula, data, ...),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(c(estimate = NA_real_, std_error = NA_real_,
               t_value = NA_real_, df = NA_real_, p_value = NA_real_,
               rho = NA_real_, pearson_ratio = NA_real_, scale = NA_real_,
               dispersion_boundary = NA_real_,
               mean_cpm = mean(counts[i, ] / totals * 1e6),
               converged = 0))
    }
    tab <- fit$coefficient_table[term, ]
    c(
      estimate = if (fit$converged) tab[["estimate"]] else NA_real_,
      std_error = tab[["std_error"]],
      t_value = tab[["t_value"]],
      df = tab[["df"]],
      p_value = tab[["p_value"]],
      rho = fit$rho,
      pearson_ratio = fit$pearson / fit$df.residual,
      scale = fit$scale,
      dispersion_boundary = as.numeric(fit$dispersion_boundary),
      mean_cpm = mean(counts[i, ] / totals * 1e6),
      converged = as.numeric(fit$converged)
    )
  }
  if (ncores > 1L && .Platform$OS.type == "unix") {
    pieces <- parallel::mclapply(
      seq_len(nrow(counts)), one_guide, mc.cores = ncores
    )
    statistics <- do.call(rbind, pieces)
  } else {
    statistics <- t(vapply(
      seq_len(nrow(counts)), one_guide, numeric(11L)
    ))
  }
  result <- data.frame(
    guide = guide,
    statistics,
    row.names = NULL,
    check.names = FALSE
  )
  result$dispersion_boundary <- as.logical(result$dispersion_boundary)
  result$converged <- as.logical(result$converged)
  if (!is.null(gene)) {
    result <- cbind(gene = gene, result)
  }
  result$fdr <- p.adjust(result$p_value, method = "BH")
  result
}

#' Calibrate guide-level t tests with negative-control guides
#'
#' Estimates a one-parameter empirical-null scale from the absolute
#' negative-control t statistics.  The scale matches their empirical
#' `(1 - alpha)` quantile to the corresponding two-sided Student t cutoff.
#' Scales below `min_scale` are truncated, so calibration need not make an
#' already conservative analysis more liberal.
#'
#' @param result A guide-level data frame returned by [bb_screen()].
#' @param control Logical vector identifying negative-control guides.
#' @param alpha Tail probability at which to estimate the null scale.
#' @param min_controls Minimum number of finite negative controls required.
#' @param min_scale Lower bound for the estimated scale.
#' @return `result` with recalibrated standard errors, t statistics, p-values,
#'   and FDR. Original inferential columns are retained with a `raw_` prefix;
#'   the scale and alpha are stored as attributes.
#' @export
bb_calibrate_controls <- function(result, control, alpha = 0.05,
                                  min_controls = 20L, min_scale = 1) {
  required <- c("estimate", "std_error", "t_value", "df", "p_value", "fdr")
  if (!is.data.frame(result) || !all(required %in% names(result))) {
    .bb_stop(
      "`result` must be a guide-level result from `bb_screen()`."
    )
  }
  if (!is.logical(control) || length(control) != nrow(result) ||
      anyNA(control)) {
    .bb_stop("`control` must be a non-missing logical vector, one per guide.")
  }
  if (length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 0.5) {
    .bb_stop("`alpha` must be one finite number between 0 and 0.5.")
  }
  if (length(min_controls) != 1L || !is.finite(min_controls) ||
      min_controls < 2) {
    .bb_stop("`min_controls` must be at least two.")
  }
  if (length(min_scale) != 1L || !is.finite(min_scale) ||
      min_scale <= 0) {
    .bb_stop("`min_scale` must be positive.")
  }

  valid <- control & is.finite(result$t_value) & is.finite(result$df)
  if (sum(valid) < as.integer(min_controls)) {
    .bb_stop(sprintf(
      "At least %d finite negative-control statistics are required.",
      as.integer(min_controls)
    ))
  }
  control_df <- unique(result$df[valid])
  if (length(control_df) != 1L || control_df <= 0) {
    .bb_stop("Finite negative-control guides must share one positive `df`.")
  }
  empirical_cutoff <- as.numeric(stats::quantile(
    abs(result$t_value[valid]),
    probs = 1 - alpha,
    names = FALSE,
    type = 8
  ))
  reference_cutoff <- stats::qt(1 - alpha / 2, df = control_df)
  scale <- max(min_scale, empirical_cutoff / reference_cutoff)

  result$raw_std_error <- result$std_error
  result$raw_t_value <- result$t_value
  result$raw_p_value <- result$p_value
  result$raw_fdr <- result$fdr
  result$std_error <- result$std_error * scale
  result$t_value <- result$t_value / scale
  result$p_value <- 2 * pt(-abs(result$t_value), df = result$df)
  result$fdr <- p.adjust(result$p_value, method = "BH")
  attr(result, "control_scale") <- scale
  attr(result, "control_alpha") <- alpha
  result
}

#' Test a shared guide effect against an empirical gene-level null
#'
#' This function is intended for exploratory screens in which several
#' independently designed guides target each gene but biological replication
#' is too limited for reliable guide-level reference distributions. It does
#' not treat guides as biological replicates. Instead, it estimates one shared
#' gene effect by inverse-variance weighting of guide coefficients, forms its
#' model-based Wald statistic, then calibrates the gene-statistic distribution
#' with a robust empirical null. It does not combine guide p-values by Fisher's
#' or Stouffer's method.
#'
#' For gene \eqn{g}, let
#' \eqn{w_{gj}=\mathrm{SE}(\widehat\beta_{gj})^{-2}}.
#' The shared effect and raw statistic are
#' \deqn{\widehat\beta_g =
#' \frac{\sum_j w_{gj}\widehat\beta_{gj}}{\sum_j w_{gj}},\qquad
#' T_g = \widehat\beta_g\sqrt{\sum_jw_{gj}}.}
#' Its null center is the median score among control genes when enough are
#' supplied, and otherwise the median across all genes. The null scale is the
#' largest of `min_scale`, the all-gene MAD, and the control-gene tail scale.
#' The all-gene MAD assumes that fewer than half of genes are active.
#'
#' The inverse-variance standard error does not model dependence among guide
#' coefficients induced by their shared sample libraries. Consequently,
#' `p_value` and `fdr` are empirical-null working quantities for ranking and
#' sensitivity analysis, not confirmatory error guarantees. Independently
#' designed guides support perturbation consistency but do not replace
#' independent biological samples.
#'
#' @param result A guide-level result from [bb_screen()]. If control
#'   calibration has already been applied, the retained `raw_std_error`
#'   column is used automatically. Nonconverged guides are excluded when a
#'   `converged` column is present.
#' @param control Optional non-missing logical vector identifying control
#'   guides. A control gene must contain only control guides.
#' @param min_guides Minimum number of finite guide scores required per gene.
#' @param alpha Tail probability used to estimate the control-gene scale.
#' @param min_control_genes Minimum number of valid control genes needed to
#'   use their median and tail scale.
#' @param min_scale Lower bound for the empirical-null scale.
#' @return A data frame with one row per testable gene. `statistic` is the
#'   empirical-null standardized gene score and `p_value` uses a standard
#'   normal working reference. Null parameters are also stored as attributes.
#' @export
bb_gene_consistency <- function(result, control = NULL, min_guides = 3L,
                                alpha = 0.05, min_control_genes = 10L,
                                min_scale = 1) {
  required <- c("gene", "estimate", "std_error")
  if (!is.data.frame(result) || !all(required %in% names(result))) {
    .bb_stop(
      "`result` must contain guide-level `gene`, `estimate`, and `std_error` columns."
    )
  }
  if (anyNA(result$gene)) {
    .bb_stop("`result$gene` cannot contain missing values.")
  }
  if (length(min_guides) != 1L || !is.finite(min_guides) ||
      min_guides < 2) {
    .bb_stop("`min_guides` must be an integer of at least two.")
  }
  if (length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 0.5) {
    .bb_stop("`alpha` must be one finite number between 0 and 0.5.")
  }
  if (length(min_control_genes) != 1L ||
      !is.finite(min_control_genes) || min_control_genes < 2) {
    .bb_stop("`min_control_genes` must be an integer of at least two.")
  }
  if (length(min_scale) != 1L || !is.finite(min_scale) ||
      min_scale <= 0) {
    .bb_stop("`min_scale` must be positive.")
  }
  min_guides <- as.integer(min_guides)
  min_control_genes <- as.integer(min_control_genes)

  if (is.null(control)) {
    control <- rep(FALSE, nrow(result))
  } else if (!is.logical(control) || length(control) != nrow(result) ||
             anyNA(control)) {
    .bb_stop("`control` must be a non-missing logical vector, one per guide.")
  }
  if ("converged" %in% names(result) &&
      (!is.logical(result$converged) || anyNA(result$converged))) {
    .bb_stop("`result$converged` must be a non-missing logical vector.")
  }
  control_by_gene <- split(control, result$gene)
  mixed_control <- vapply(
    control_by_gene,
    function(value) any(value) && !all(value),
    logical(1L)
  )
  if (any(mixed_control)) {
    .bb_stop("A gene cannot mix control and non-control guides.")
  }

  standard_error <- if ("raw_std_error" %in% names(result)) {
    result$raw_std_error
  } else {
    result$std_error
  }
  valid <- is.finite(result$estimate) &
    is.finite(standard_error) &
    standard_error > 0
  if ("converged" %in% names(result)) {
    valid <- valid & result$converged
  }
  groups <- split(seq_len(nrow(result)), result$gene)
  pieces <- lapply(names(groups), function(gene_name) {
    all_index <- groups[[gene_name]]
    index <- all_index[valid[all_index]]
    if (length(index) < min_guides) {
      return(NULL)
    }
    guide_weight <- 1 / standard_error[index]^2
    gene_estimate <- sum(
      guide_weight * result$estimate[index]
    ) / sum(guide_weight)
    gene_standard_error <- sqrt(1 / sum(guide_weight))
    nonzero <- result$estimate[index] != 0
    agreement <- if (gene_estimate == 0 || !any(nonzero)) {
      NA_real_
    } else {
      mean(
        sign(result$estimate[index][nonzero]) == sign(gene_estimate)
      )
    }
    data.frame(
      gene = gene_name,
      n_guides = length(index),
      estimate = gene_estimate,
      std_error = gene_standard_error,
      raw_statistic = gene_estimate / gene_standard_error,
      guide_direction_agreement = agreement,
      converged_fraction = if ("converged" %in% names(result)) {
        mean(result$converged[all_index])
      } else {
        NA_real_
      },
      control_gene = all(control[all_index]),
      row.names = NULL
    )
  })
  pieces <- Filter(Negate(is.null), pieces)
  if (length(pieces) < 2L) {
    .bb_stop("At least two genes must have enough finite guide scores.")
  }
  gene_result <- do.call(rbind, pieces)
  rownames(gene_result) <- NULL

  global_center <- stats::median(gene_result$raw_statistic)
  global_scale <- stats::mad(
    gene_result$raw_statistic,
    center = global_center,
    constant = 1 / stats::qnorm(0.75)
  )
  if (!is.finite(global_scale) || global_scale <= 0) {
    global_scale <- 1
  }

  control_statistic <- gene_result$raw_statistic[
    gene_result$control_gene
  ]
  enough_controls <- length(control_statistic) >= min_control_genes
  null_center <- if (enough_controls) {
    stats::median(control_statistic)
  } else {
    global_center
  }
  control_scale <- if (enough_controls) {
    as.numeric(stats::quantile(
      abs(control_statistic - null_center),
      probs = 1 - alpha,
      names = FALSE,
      type = 8
    )) / stats::qnorm(1 - alpha / 2)
  } else {
    NA_real_
  }
  scale_candidates <- c(min_scale, global_scale, control_scale)
  null_scale <- max(scale_candidates[is.finite(scale_candidates)])

  gene_result$statistic <-
    (gene_result$raw_statistic - null_center) / null_scale
  gene_result$p_value <- 2 * stats::pnorm(-abs(gene_result$statistic))
  gene_result$fdr <- stats::p.adjust(gene_result$p_value, method = "BH")
  attr(gene_result, "null_center") <- null_center
  attr(gene_result, "null_scale") <- null_scale
  attr(gene_result, "global_scale") <- global_scale
  attr(gene_result, "control_scale") <- control_scale
  attr(gene_result, "control_genes") <- length(control_statistic)
  attr(gene_result, "null_assumption") <-
    paste(
      "Shared-effect guide-consistency empirical null;",
      "not biological-replicate inference."
    )
  gene_result
}
