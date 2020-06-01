#' Doubly robust estimates of for evaluating effects 
#' of treatments on ordinal outcomes.
#' 
#' The available parameters for evaluating treatment efficacy are:
#' \itemize{
#' \item Difference in (weighted) means: The outcome levels are treated
#' numerically, with each level possibly assigned a weight. The difference
#' in average outcomes is computed.
#' \item Log odds ratio: The comparison describes the average log-odds (treatment
#' level 1 versus 0) of the cumulative probability for each level of the outcome.
#' \item Mann-Whitney: The probability that a randomly-selected individual 
#' receiving treatment 1 will have a larger outcome value than a randomly selected
#' individual receiving treatment 0 (with ties assigned weight 1/2). 
#' }
#' 
#' In each case, estimates are constructed by obtaining a doubly robust 
#' estimate of the cumulative distribution function (CDF) for each treatment group. 
#' This is achieved by fitting a (working) proportional odds model that includes
#' inverse probability of treatment weights. The inclusion of these weights 
#' ensures that, so long as the working model includes intercept terms,
#' the resultant estimate of the CDF is an augmented inverse 
#' probability of treatment weighted estimate. This implies that the estimate is
#' nonparametric efficient if the working model contains the truth; however, 
#' even if the working model does not contain the truth, the CDF estimates are 
#' consistent and asymptotically normal with variance expected to dominate that
#' of an unadjusted estimate of the same treatment effect. 
#' 
#' The CDF estimates are subsequently mapped into estimates of each requested 
#' parameter for evaluating treatment effects. The double robustness and efficiency 
#' properties of the CDF estimates extend to these quantities as well. Confidence
#' intervals and hypothesis tests can be carried out in closed form using Wald-style
#' intervals and tests or using a nonparametric corrected and accelerated bootstrap 
#' (BCa). Inference for the CDF and probability mass function is also returned and 
#' can be used for subsequent visualizations (see \code{plot.drord}). 
#' 
#' @param out A \code{numeric} vector containing the outcomes. Missing outcomes are 
#' allowed. 
#' @param treat A \code{numeric} vector containing treatment status. Missing
#' values are not allowed unless the corresponding entry in \code{out} is also missing. 
#' Only values of 0 or 1 are treated as actual treatment levels. Any other value is assumed 
#' to encode a value for which the outcome is missing and the corresponding outcome value is 
#' ignored. 
#' @param covar A \code{data.frame} containing the covariates to include in the working
#' proportional odds model. 
#' @param out_levels A \code{numeric} vector containing all ordered levels of the 
#' outcome. 
#' @param out_form The right-hand side of a regression formula for the working proportional 
#' odds model. NOTE: THIS FORMULA MUST NOT SUPPRESS THE INTERCEPT. 
#' @param out_model Which R function should be used to fit the proportional odds 
#' model. Options are \code{"polr"} (from the \code{MASS} package), 
#' "vglm" (from the \code{VGAM} package), or \code{"clm"} (from the \code{ordinal} package).
#' @param out_weights A vector of \code{numeric} weights with length equal to the length 
#' of \code{out_levels}. 
#' @param treat_form The right-hand side of a regression formula for the working model of
#' treatment probability as a function of covariates
#' @param param A vector of \code{characters} indicating which of the three treatment
#' effect parameters should be estimated (\code{"weighted_mean", "log_odds",} 
#' and/or \code{"mann_whitney"}).
#' @param ci A vector of \code{characters} indicating which confidence intervals
#' should be computed (\code{"bca"} and/or \code{"wald"}) 
#' @param alpha Confidence intervals have nominal level 1-\code{alpha}. 
#' @param nboot Number of bootstrap replicates used to compute bootstrap confidence
#' intervals. 
#' @param return_models If \code{TRUE} the fitted working proportional odds models
#' and treatment probability models are returned. 
#' @param stratify If \code{TRUE}, then a fully stratified estimator is computed, i.e.,
#' the empirical CDF of each treatment arm is estimated stratifying by levels of 
#' \code{covar}. For now, this option is limited to univariate covariates. 
#' @param est_dist A \code{boolean} indicating whether estimates of the CDF and PMF 
#' should be computed and returned. For real data analysis, we generally recommend 
#' leaving as \code{TRUE}; however, when studying performance in simulations, it can 
#' save time to set to \code{FALSE}. 
#' @param ... Other options (not currently used). 
#' @export
#' 
#' @return An object of class \code{drord}. In addition to information related to 
#' how \code{drord} was called, the output contains the following:
#' \describe{
#'   \item{log_odds}{inference pertaining to the log-odds parameter. \code{NULL} if 
#' 					 this parameter not requested in call to \code{drord}.}
#'   \item{mann_whitney}{inference pertaining to the Mann-Whitney parameter. \code{NULL} if 
#' 					 this parameter not requested in call to \code{drord}.}
#'   \item{weighted_mean}{inference pertaining to weighted mean parameter. \code{NULL} if 
#' 					 this parameter not requested in call to \code{drord}.}
#'   \item{cdf}{inference pertaining to the treatment-specific CDFs. See the \code{plot}
#' 				method for a convenient way of visualizing this information. \code{NULL} 
#' 				if \code{est_dist = FALSE} in call to \code{drord}.}
#' 	 \item{pmf}{inference pertaining to the treatment-specific PMFs. See the \code{plot}
#' 				method for a convenient way of visualizing this information. \code{NULL} 
#' 				if \code{est_dist = FALSE} in call to \code{drord}.}
#'   \item{treat_mod}{the fitted model for the probability of treatment as a function
#' 					  of covariates. \code{NULL} if \code{return_models = FALSE}}
#'   \item{out_mod}{the proportional odds model fit in each treatment arm. named entries 
#' 					in list indicate the corresponding treatment arm. 
#' 					\code{NULL} if \code{return_models = FALSE} or \code{stratify = TRUE}.}
#' }
#'@examples 
#' data(covid19)
#' 
#' # get estimates of all parameters based on main-effects 
#' # proportional odds model and intercept-only propensity model
#' fit <- drord(out = covid19$out, treat = covid19$treat, 
#'              covar = covid19[, "age_grp", drop = FALSE])
#' 
#' # get estimates of all parameters based on proportional odds and
#' # propensity model that treats age_grp as categorical
#' fit2 <- drord(out = covid19$out, treat = covid19$treat, 
#'               covar = covid19[, "age_grp", drop = FALSE],
#' 				 out_form = "factor(age_grp)",
#' 				 treat_form = "factor(age_grp)")
#' 
#' # obtain estimator stratified by age group 
#' fit3 <- drord(out = covid19$out, treat = covid19$treat, 
#'               covar = covid19[, "age_grp", drop = FALSE],
#' 				 stratify = TRUE)
#' 
#' # demonstration with missing outcome data
#' covid19$out[1:5] <- NA
#' 
#' # propensity model should now adjust for covariates to address
#' # the potential for informative missingness
#' fit4 <- drord(out = covid19$out, treat = covid19$treat, 
#'               covar = covid19[, "age_grp", drop = FALSE],
#' 				 treat_form = "age_grp")

drord <- function(
  out,
  treat,
  covar,
  out_levels = sort(unique(out)), # must be in order for this to work!
  out_form = ".",
  out_weights = rep(1, length(out_levels)),
  out_model = c("clm", "polr", "vglm"),
  treat_form = "1",
  param = c("weighted_mean", "log_odds", "mann_whitney"),
  ci = "wald",
  alpha = 0.05, 
  nboot = 1e4,
  return_models = TRUE, 
  est_dist = TRUE, 
  stratify = FALSE, 
  ...
){

	# recode treat for NA outcomes
	if(any(is.na(out))){
		treat[is.na(out)] <- -1
	}

	# obtain estimate of treatment probabilities
	treat_prob_fit <- estimate_treat_prob(treat = treat,
	                                      covar = covar,
	                                      treat_form = treat_form,
	                                      return_models = TRUE)
	treat_prob_est <- treat_prob_fit$gn
	if(return_models){
		treat_mod <- treat_prob_fit$fm
	}else{
		treat_mod <- NULL
	}

	# obtain estimate of conditional PMF for each treatment level
	pmf_fit <- estimate_pmf(out = out, treat = treat, 
	                        covar = covar, out_levels = out_levels,
	                        out_form = out_form, treat_prob_est = treat_prob_est,
	                        out_model = out_model, 
	                        stratify = stratify, return_models = return_models)
  	pmf_est <- pmf_fit$pmf
  	if(return_models){
  		out_mod <- pmf_fit$fm
  	}else{
  		out_mod <- NULL
  	}

  	# map into estimates of CDF
	cdf_est <- estimate_cdf(pmf_est = pmf_est)

	# get estimate of weighted mean
	if("weighted_mean" %in% param){
		wmean_est <- estimate_wmean(
          pmf_est = pmf_est, treat = treat, out = out, out_levels = out_levels, 
          out_weights = out_weights, treat_prob_est = treat_prob_est, 
          return_cov = "wald" %in% ci
        )
        wmean_ci <- estimate_ci_wmean(out = out, treat = treat,
                                      covar = covar, wmean_est = wmean_est,
                                      alpha = alpha, out_levels = out_levels,
                                      out_form = out_form, out_weights = out_weights,
                                      out_model = out_model,
                                      treat_form = treat_form, ci = ci, nboot = nboot)
	}else{
		wmean_est <- NULL
		wmean_ci <- NULL
	}

	# get estimate of log_odds ratio
	if("log_odds" %in% param){
		logodds_est <- estimate_logodds(cdf_est = cdf_est)
		logodds_ci <- estimate_ci_logodds(logodds_est = logodds_est, 
		                                  cdf_est = cdf_est,
		                                  treat_prob_est = treat_prob_est, 
		                                  treat_form = treat_form,
		                                  out_form = out_form,
		                                  out_model = out_model, 
		                                  treat = treat, ci = ci, 
		                                  out = out, alpha = alpha,
		                                  nboot = nboot, covar = covar, 
		                                  out_levels = out_levels)
	}else{
		logodds_est <- NULL
		logodds_ci <- NULL
	}

	if("mann_whitney" %in% param){
		mannwhitney_est <- estimate_mannwhitney(cdf_est = cdf_est,
		                                        pmf_est = pmf_est)
		mannwhitney_ci <- estimate_ci_mannwhitney(mannwhitney_est = mannwhitney_est, 
                                  cdf_est = cdf_est,
                                  pmf_est = pmf_est,
                                  treat_prob_est = treat_prob_est, 
                                  treat_form = treat_form, 
                                  out_form = out_form,
                                  out_model = out_model,
                                  treat = treat, ci = ci, 
                                  out = out, alpha = alpha,
                                  nboot = nboot, out_levels = out_levels,
                                  covar = covar)
	}else{
		mannwhitney_est <- NULL
		mannwhitney_ci <- NULL
	}

	if(est_dist){
		# get pt-wise and simultaneous confidence intervals for cdf
		marg_cdf_est <- marginalize_cdf(cdf_est = cdf_est)
		marg_pmf_est <- marginalize_pmf(pmf_est = pmf_est)

		marg_dist_ci <- estimate_ci_marg_dist(marg_cdf_est = marg_cdf_est,
		                                      marg_pmf_est = marg_pmf_est,
		                                      out_levels = out_levels,
		                                      cdf_est = cdf_est, pmf_est = pmf_est,
		                             		  treat_prob_est = treat_prob_est, 
		                             		  treat_form = treat_form,
		                             		  out_form = out_form, 
		                             		  out_model = out_model,
		                             		  treat = treat, ci = ci, 
		                             		  out = out, alpha = alpha,
		                             		  covar = covar, 
		                             		  nboot = nboot)
	}else{
		marg_cdf_est <- NULL; marg_pmf_est <- NULL
		marg_dist_ci <- list(cdf = NULL, pmf = NULL)
	}

	# structure output
	rout <- list(log_odds = list(est = logodds_est, ci = logodds_ci),
	            mann_whitney = list(est = mannwhitney_est, ci = mannwhitney_ci),
	            weighted_mean = list(est = wmean_est, ci = wmean_ci),
	            cdf = list(est = marg_cdf_est, ci = marg_dist_ci$cdf),
	            pmf = list(est = marg_pmf_est, ci = marg_dist_ci$pmf),
	            out_levels = out_levels, treat_mod = treat_mod,
	            out_mod = out_mod, param = param, 
	            ci = ci, alpha = alpha, est_dist = est_dist)
	class(rout) <- "drord"
	return(rout)
}

