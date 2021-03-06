% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dist_fn.R
\name{evaluate_marg_pmf_eif}
\alias{evaluate_marg_pmf_eif}
\title{Get eif estimates for treatment-specific PMF}
\usage{
evaluate_marg_pmf_eif(pmf_est, treat_prob_est, treat, out, out_levels)
}
\arguments{
\item{pmf_est}{Estimated conditional PMF for \code{trt_level}.}

\item{treat_prob_est}{Estimated propensity for \code{trt_level}.}

\item{treat}{A \code{numeric} vector containing treatment status. Should only assume 
a value 0 or 1.}

\item{out}{A \code{numeric} vector containing the outcomes. Missing outcomes are 
allowed.}

\item{out_levels}{A \code{numeric} vector containing all ordered levels of the 
outcome.}
}
\value{
a list of eif estimates
}
\description{
Get eif estimates for treatment-specific PMF
}
