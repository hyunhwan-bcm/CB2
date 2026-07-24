#!/usr/bin/env Rscript

# Representative performance check for the compiled beta-binomial IRLS kernel
# and an end-to-end guide screen. Run after installing CB2:
#
#   Rscript inst/benchmarks/bbreg.R

library(CB2)

set.seed(20260723)
n_guides <- 10000L
design <- data.frame(
  dose = rep(seq(-1, 1, length.out = 8), 2),
  batch = factor(rep(c("a", "b"), each = 8))
)
totals <- rep(60000, nrow(design))
eta <- outer(runif(n_guides, -7.2, -6.6), rep(1, nrow(design))) +
  outer(runif(n_guides, -0.4, 0.4), design$dose)
mu <- plogis(eta)
rho <- 0.0015
precision <- 1 / rho - 1
latent <- matrix(
  rbeta(length(mu), mu * precision, (1 - mu) * precision),
  nrow = n_guides
)
counts <- matrix(
  rbinom(length(mu), rep(totals, each = n_guides), latent),
  nrow = n_guides
)

elapsed <- system.time({
  result <- bb_screen(
    counts, design, ~ dose + batch, "dose", totals = totals,
    ncores = min(4L, parallel::detectCores(logical = FALSE))
  )
})[["elapsed"]]

cat(sprintf(
  paste0(
    "%d guides in %.3f seconds (%.1f guides/second); ",
    "%.1f%% converged.\n"
  ),
  n_guides, elapsed, n_guides / elapsed, 100 * mean(result$converged)
))
