\name{bb_calibrate_controls}
\alias{bb_calibrate_controls}
\title{Calibrate guide-level t tests with negative-control guides}
\usage{
bb_calibrate_controls(
  result,
  control,
  alpha = 0.05,
  min_controls = 20L,
  min_scale = 1
)
}
\arguments{
\item{result}{A guide-level data frame returned by \code{bb_screen()}.}

\item{control}{Logical vector identifying negative-control guides.}

\item{alpha}{Tail probability at which to estimate the null scale.}

\item{min_controls}{Minimum number of finite negative controls required.}

\item{min_scale}{Lower bound for the estimated scale.}
}
\value{
\code{result} with recalibrated standard errors, t statistics, p-values, and
FDR. Original inferential columns are retained with a \code{raw_} prefix; the
scale and alpha are stored as attributes.
}
\description{
Estimates a one-parameter empirical-null scale from the absolute
negative-control t statistics. The scale matches their empirical
\code{(1 - alpha)} quantile to the corresponding two-sided Student t cutoff.
Scales below \code{min_scale} are truncated, so calibration need not make an
already conservative analysis more liberal.
}
