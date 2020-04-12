#' Obtain estimates of ...
#' 
#' @param ...
#' @export
#'@examples 
#' n <- 200
#' covar <- data.frame(x1 = runif(n), x2 = runif(n))
#' treat <- rbinom(n, 1, 1/2)
#' out <- rbinom(n, 3, plogis(-1 + covar$x1 - covar$x2 + 0.2 * treat))
#' out_levels <- 0:3
#' out_form <- "x1 + x2"
#' out_weights <- rep(1, 4)
#' treat_form <- "x1 + x2"
#' ci <- c("bca", "wald")
#' test <- c("perm", "wald")
#' alpha <- 0.05 
#' nboot <- 1e2
#' out_levels <- 0:3

drord <- function(
  out,
  treat,
  covar,
  out_levels = order(unique(out)), # must be in order for this to work!
  out_form = NULL,
  out_weights = rep(1, length(out_levels)),
  treat_form = "1",
  param = c("weighted-mean", "log-odds", "mann-whitney"),
  ci = c("bca", "wald"),
  alpha = 0.05, 
  test = c("perm", "wald"), # not implemented yet
  nboot = 1e4,
  nperm = 1e4, # not implemented yet
  est_dist = FALSE, # temporary to make sims go faster under time crunch
  ...
){
	# obtain estimate of treatment probabilities
	treat_prob_est <- estimate_treat_prob(treat = treat,
	                                      covar = covar,
	                                      treat_form = treat_form)

	# obtain estimate of conditional PMF for each treatment level
	pmf_est <- estimate_pmf(out = out, treat = treat, 
	                        covar = covar, out_levels = out_levels,
	                        out_form = out_form, treat_prob_est = treat_prob_est)
  	
  	# map into estimates of CDF
	cdf_est <- estimate_cdf(pmf_est = pmf_est)

	# get estimate of weighted mean
	if("weighted-mean" %in% param){
		wmean_est <- estimate_wmean(
          pmf_est = pmf_est, treat = treat, out = out, out_levels = out_levels, 
          out_weights = out_weights, treat_prob_est = treat_prob_est, 
          return_cov = "wald" %in% ci | "wald" %in% test
        )
        wmean_ci <- estimate_ci_wmean(out = out, treat = treat,
                                      covar = covar, wmean_est = wmean_est,
                                      alpha = alpha, out_levels = out_levels,
                                      out_form = out_form, out_weights = out_weights,
                                      treat_form = treat_form, ci = ci, nboot = nboot)
	}else{
		wmean_est <- NULL
		wmean_ci <- NULL
	}

	# get estimate of log-odds ratio
	if("log-odds" %in% param){
		logodds_est <- estimate_logodds(cdf_est = cdf_est)
		logodds_ci <- estimate_ci_logodds(logodds_est = logodds_est, 
		                                  cdf_est = cdf_est,
		                                  treat_prob_est = treat_prob_est, 
		                                  treat_form = treat_form,
		                                  out_form = out_form,
		                                  treat = treat, ci = ci, 
		                                  out = out, alpha = alpha,
		                                  nboot = nboot, covar = covar, 
		                                  out_levels = out_levels)
	}else{
		logodds_est <- NULL
		logodds_ci <- NULL
	}

	if("mann-whitney" %in% param){
		mannwhitney_est <- estimate_mannwhitney(cdf_est = cdf_est,
		                                        pmf_est = pmf_est)
		mannwhitney_ci <- estimate_ci_mannwhitney(mannwhitney_est = mannwhitney_est, 
                                  cdf_est = cdf_est,
                                  pmf_est = pmf_est,
                                  treat_prob_est = treat_prob_est, 
                                  treat_form = treat_form, 
                                  out_form = out_form,
                                  treat = treat, ci = ci, 
                                  out = out, alpha = alpha,
                                  nboot = nboot, out_levels = out_levels)
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
		                             		  treat = treat, ci = ci, 
		                             		  out = out, alpha = alpha,
		                             		  nboot = nboot)
	}else{
		marg_cdf_est <- NULL; marg_pmf_est <- NULL
		marg_dist_ci <- list(cdf = NULL, pmf = NULL)
	}

	# structure output
	out <- list(logodds = list(est = logodds_est, ci = logodds_ci),
	            mannwhitney = list(est = mannwhitney_est, ci = mannwhitney_ci),
	            wmean = list(est = wmean_est, ci = wmean_ci),
	            cdf = list(est = marg_cdf_est, ci = marg_dist_ci$cdf),
	            pmf = list(est = marg_pmf_est, ci = marg_dist_ci$pmf))
	return(out)
}

