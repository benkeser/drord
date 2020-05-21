#' Compute confidence interval/s for the treatment specific 
#' PMF and CDF.
#' 
#' @param marg_cdf_est Point estimate of treatment-specific CDF.
#' @param marg_pmf_est Point estimate of treatment-specific PMF.
#' @param cdf_est Estimates of treatment-specific conditional CDF.
#' @param pmf_est Estimates of treatment-specific conditional PMF.
#' @param out A \code{numeric} vector containing the outcomes. Missing outcomes are 
#' allowed. 
#' @param treat A \code{numeric} vector containing treatment status. Missing
#' values are not allowed unless the corresponding entry in \code{out} is also missing. 
#' Only values of 0 or 1 are treated as actual treatment levels. Any other value is assumed 
#' to encode a value for which the outcome is missing and the corresponding outcome value is 
#' ignored. 
#' @param covar A \code{data.frame} containing the covariates to include in the working
#' proportional odds model. 
#' @param alpha Confidence intervals have nominal level 1-\code{alpha}. 
#' @param out_levels A \code{numeric} vector containing all ordered levels of the 
#' outcome. 
#' @param out_form The right-hand side of a regression formula for the working proportional 
#' odds model. NOTE: THIS FORMULA MUST NOT SUPPRESS THE INTERCEPT. 
#' @param out_model Which R function should be used to fit the proportional odds 
#' model. Options are \code{"polr"} (from the \code{MASS} package), 
#' "vglm" (from the \code{VGAM} package), or \code{"clm"} (from the \code{ordinal} package).
#' @param treat_form The right-hand side of a regression formula for the working model of
#' treatment probability as a function of covariates
#' @param ci A vector of \code{characters} indicating which confidence intervals
#' should be computed (\code{"bca"} and/or \code{"wald"}) 
#' @param nboot Number of bootstrap replicates used to compute bootstrap confidence
#' intervals. 
#' @param treat_prob_est Estimated probability of treatments, output from call
#' to \code{estimate_treat_prob}.
#' @return List of lists (\code{cdf} and \code{pmf}) with \code{wald} and \code{bca}-estimated confidence 
#' intervals for the marginal treatment-specific distribution functions. 

estimate_ci_marg_dist <- function(marg_cdf_est,
                                  marg_pmf_est,                                  
                                  cdf_est,
                                  pmf_est,
                                  covar,
	                              treat_prob_est, 
	                              treat_form, out_form,
	                              treat, ci, out_levels, 
	                              out_model,
	                              out, alpha,
	                              nboot){
	if("wald" %in% ci){
		# pointwise ci
		marg_cdf_eif <- evaluate_marg_cdf_eif(cdf_est = cdf_est, 
		                                   treat_prob_est = treat_prob_est, 
		                                   treat = treat, out = out,
		                                   out_levels = out_levels)
		marg_cdf_ptwise_wald_ci <- evaluate_marg_cdf_ptwise_ci(marg_cdf_est = marg_cdf_est,
		                                                  marg_cdf_eif = marg_cdf_eif,
		                                                  alpha = alpha)
		# simultaneous CI
		marg_cdf_simul_wald_ci <- evaluate_marg_dist_simul_ci(marg_dist_est = marg_cdf_est,
	                                                    marg_dist_eif = marg_cdf_eif,
	                                                    alpha = alpha, 
	                                                    remove_last = TRUE)

		marg_pmf_eif <- evaluate_marg_pmf_eif(pmf_est = pmf_est, 
		                                   treat_prob_est = treat_prob_est, 
		                                   treat = treat, out = out,
		                                   out_levels = out_levels)

		marg_pmf_ptwise_wald_ci <- evaluate_marg_pmf_ptwise_ci(marg_pmf_est = marg_pmf_est,
		                                                  marg_pmf_eif = marg_pmf_eif,
		                                                  alpha = alpha)
		# simultaneous CI
		marg_pmf_simul_wald_ci <- evaluate_marg_dist_simul_ci(marg_dist_est = marg_pmf_est,
	                                                    marg_dist_eif = marg_pmf_eif,
	                                                    alpha = alpha, 
	                                                    remove_last = FALSE)

		wald_ci_cdf <- list(list(ptwise = marg_cdf_ptwise_wald_ci[[1]], simul = marg_cdf_simul_wald_ci[[1]]),
		                           list(ptwise = marg_cdf_ptwise_wald_ci[[2]], simul = marg_cdf_simul_wald_ci[[2]]))
		wald_ci_pmf <- list(list(ptwise = marg_pmf_ptwise_wald_ci[[1]], simul = marg_pmf_simul_wald_ci[[1]]),
		                           list(ptwise = marg_pmf_ptwise_wald_ci[[2]], simul = marg_pmf_simul_wald_ci[[2]]))
	}else{
		wald_ci_cdf <- list(ptwise = NULL,
		                	simul = NULL)
		wald_ci_pmf <- list(ptwise = NULL,
		                	simul = NULL)
		
	}

	if("bca" %in% ci){
		bca_ci <- bca_marg_dist(treat = treat, covar = covar,
		                       out = out, nboot = nboot, 
		                       treat_form = treat_form, 
		                       out_levels = out_levels, 
		                       out_form = out_form,
		                       marg_cdf_est = marg_cdf_est,
		                       marg_pmf_est = marg_pmf_est,
		                       alpha = alpha,
		                       out_model = out_model)
		bca_ci_cdf <- bca_ci$cdf
		bca_ci_pmf <- bca_ci$pmf
	}else{
		bca_ci_cdf <- NULL
		bca_ci_pmf <- NULL
	}
	return(list(cdf = list(wald = wald_ci_cdf, bca = bca_ci_cdf),
	       		pmf = list(wald = wald_ci_pmf, bca = bca_ci_pmf)))
}

#' Marginalize over empirical distribution to obtain marginal
#' treatment-specific CDF estimate.
#' 
#' @param cdf_est Estimates of treatment-specific conditional CDF.

marginalize_cdf <- function(cdf_est){
	lapply(cdf_est, colMeans)
}



#' Compute a BCa bootstrap confidence interval for the weighted mean. The code is 
#' based on the slides found here: http://users.stat.umn.edu/~helwig/notes/bootci-Notes.pdf
#' 
#' @param out A \code{numeric} vector containing the outcomes. Missing outcomes are 
#' allowed. 
#' @param treat A \code{numeric} vector containing treatment status. Should only assume 
#' a value 0 or 1. 
#' @param covar A \code{data.frame} containing the covariates to include in the working
#' proportional odds model. 
#' @param nboot Number of bootstrap replicates used to compute bootstrap confidence
#' intervals. 
#' @param treat_form The right-hand side of a regression formula for the working model of
#' treatment probability as a function of covariates
#' @param out_levels A \code{numeric} vector containing all ordered levels of the 
#' outcome. 
#' @param out_form The right-hand side of a regression formula for the working proportional 
#' odds model. NOTE: THIS FORMULA MUST NOT SUPPRESS THE INTERCEPT. 
#' @param out_model Which R function should be used to fit the proportional odds 
#' model. Options are \code{"polr"} (from the \code{MASS} package), 
#' "vglm" (from the \code{VGAM} package), or \code{"clm"} (from the \code{ordinal} package).
#' @param marg_cdf_est Point estimate of treatment-specific CDF.
#' @param marg_pmf_est Point estimate of treatment-specific PMF.
#' @param alpha Level of confidence interval.
#' @return List (\code{cdf}, \code{pmf}) of lists (\code{treat=1}, \code{treat=0}) of
#' confidence intervals for distributions.
bca_marg_dist <- function(treat, covar, out, nboot, 
                      treat_form, out_levels, out_form, out_model,
                      marg_cdf_est, marg_pmf_est, alpha = 0.05){
	K <- length(out_levels)
	boot_samples <- replicate(nboot, 
	                          one_boot_marg_dist(treat = treat, 
                                               covar = covar, 
	                                           out = out, 
	                                           treat_form = treat_form, 
	                                           out_levels = out_levels, 
	                                           out_form = out_form,
	                                           out_model = out_model))

	jack_samples <- jack_marg_cdf(treat = treat,
	                           covar = covar,
	                           out = out, 
	                           treat_form = treat_form, 
	                           out_levels = out_levels, 
	                           out_form = out_form,
	                           out_model = out_model)

	ci_cdf_1 <- compute_trt_spec_bca_intervals(dist = "cdf",
	                                           trt = 1, 
	                                           marg_est = marg_cdf_est,
	                                           boot_samples = boot_samples,
	                                           jack_samples = jack_samples,
	                                           alpha = alpha)
	ci_cdf_0 <- compute_trt_spec_bca_intervals(dist = "cdf",
	                                           trt = 0, 
	                                           marg_est = marg_cdf_est,
	                                           boot_samples = boot_samples,
	                                           jack_samples = jack_samples,
	                                           alpha = alpha)
	ci_pmf_1 <- compute_trt_spec_bca_intervals(dist = "pmf",
	                                           trt = 1, 
	                                           marg_est = marg_pmf_est,
	                                           boot_samples = boot_samples,
	                                           jack_samples = jack_samples,
	                                           alpha = alpha)
	ci_pmf_0 <- compute_trt_spec_bca_intervals(dist = "pmf",
	                                           trt = 0, 
	                                           marg_est = marg_pmf_est,
	                                           boot_samples = boot_samples,
	                                           jack_samples = jack_samples,
	                                           alpha = alpha)
	
	rslt <- list(cdf = list(ci_cdf_1, ci_cdf_0),
	             pmf = list(ci_pmf_1, ci_pmf_0))
	return(rslt)
}

#' Compute jackknife distribution estimates.
#' @param out A \code{numeric} vector containing the outcomes. Missing outcomes are 
#' allowed. 
#' @param treat A \code{numeric} vector containing treatment status. Should only assume 
#' a value 0 or 1. 
#' @param covar A \code{data.frame} containing the covariates to include in the working
#' proportional odds model. 
#' @param treat_form The right-hand side of a regression formula for the working model of
#' treatment probability as a function of covariates
#' @param out_levels A \code{numeric} vector containing all ordered levels of the 
#' outcome. 
#' @param out_form The right-hand side of a regression formula for the working proportional 
#' odds model. NOTE: THIS FORMULA MUST NOT SUPPRESS THE INTERCEPT. 
#' @param out_model Which R function should be used to fit the proportional odds 
#' model. Options are \code{"polr"} (from the \code{MASS} package), 
#' "vglm" (from the \code{VGAM} package), or \code{"clm"} (from the \code{ordinal} package).
#' @return Jackknife estimated distributions
#' 
jack_marg_cdf <- function(treat, covar, out, treat_form, 
                          out_levels, out_form, out_model){
  	marg_cdf_jack_est <- sapply(seq_along(out), function(i){
		marg_cdf_minusi <- get_one_marg_dist(treat = treat[-i],
		                                covar = covar[-i, , drop = FALSE],
		                                out = out[-i],
		                                treat_form = treat_form,
		                                out_levels = out_levels,
		                                out_form = out_form,
		                                out_model = out_model)  		
		return(marg_cdf_minusi)
  	})
  	return(marg_cdf_jack_est)
}

#' Get one bootstrap computation of the CDF and PMF estimates 
#' 
#' @param out A \code{numeric} vector containing the outcomes. Missing outcomes are 
#' allowed. 
#' @param treat A \code{numeric} vector containing treatment status. Should only assume 
#' a value 0 or 1. 
#' @param covar A \code{data.frame} containing the covariates to include in the working
#' proportional odds model. 
#' @param treat_form The right-hand side of a regression formula for the working model of
#' treatment probability as a function of covariates
#' @param out_levels A \code{numeric} vector containing all ordered levels of the 
#' outcome. 
#' @param out_form The right-hand side of a regression formula for the working proportional 
#' odds model. NOTE: THIS FORMULA MUST NOT SUPPRESS THE INTERCEPT. 
#' @param out_model Which R function should be used to fit the proportional odds 
#' model. Options are \code{"polr"} (from the \code{MASS} package), 
#' "vglm" (from the \code{VGAM} package), or \code{"clm"} (from the \code{ordinal} package).
#' @return Estimates of CDF and PMF for a particular bootstrap sample. 

one_boot_marg_dist <- function(treat, covar, out, treat_form, 
                              out_levels, out_form, out_model){
	boot_idx <- sample(seq_along(out), replace = TRUE)
	marg_cdf_boot_est <- get_one_marg_dist(treat = treat[boot_idx],
	                                covar = covar[boot_idx, , drop = FALSE],
	                                out = out[boot_idx],
	                                treat_form = treat_form,
	                                out_levels = out_levels,
	                                out_form = out_form,
	                                out_model = out_model)
	return(marg_cdf_boot_est)
}


#' Compute one estimate of the marginal CDF/PMF on a given data set. 
#' 
#' @param out A \code{numeric} vector containing the outcomes. Missing outcomes are 
#' allowed. 
#' @param treat A \code{numeric} vector containing treatment status. Should only assume 
#' a value 0 or 1. 
#' @param covar A \code{data.frame} containing the covariates to include in the working
#' proportional odds model. 
#' @param treat_form The right-hand side of a regression formula for the working model of
#' treatment probability as a function of covariates
#' @param out_levels A \code{numeric} vector containing all ordered levels of the 
#' outcome. 
#' @param out_form The right-hand side of a regression formula for the working proportional 
#' odds model. NOTE: THIS FORMULA MUST NOT SUPPRESS THE INTERCEPT. 
#' @param out_model Which R function should be used to fit the proportional odds 
#' model. Options are \code{"polr"} (from the \code{MASS} package), 
#' "vglm" (from the \code{VGAM} package), or \code{"clm"} (from the \code{ordinal} package).
#' @return List of estimated cdf/pmf for these input data. 
get_one_marg_dist <- function(treat, covar, treat_form, out_model,
                             out, out_levels, out_form){
	# obtain estimate of treatment probabilities
	treat_prob_fit <- estimate_treat_prob(treat = treat,
	                                      covar = covar,
	                                      treat_form = treat_form,
	                                      return_models = FALSE)
	treat_prob_est <- treat_prob_fit$gn

	# obtain estimate of conditional PMF for each treatment level
	pmf_fit <- estimate_pmf(out = out, treat = treat, 
	                        covar = covar, out_levels = out_levels,
	                        out_form = out_form, treat_prob_est = treat_prob_est,
	                        out_model = out_model,
	                        return_models = FALSE)
  	pmf_est <- pmf_fit$pmf

	cdf_est <- estimate_cdf(pmf_est = pmf_est)

  	marg_cdf_est <- marginalize_cdf(cdf_est = cdf_est)
  	marg_pmf_est <- marginalize_pmf(pmf_est = pmf_est)

  	return(list(cdf = marg_cdf_est, pmf = marg_pmf_est))
}

#' Evaluate simultaneous confidence interval for marginal PMF or CDF. 
#' @param marg_dist_est The point estimate of the marginal CDF/PMF distribution
#' @param marg_dist_eif The EIF estimates for the marginal CDF/PMF estimates
#' @param alpha Confidence intervals have nominal level 1-\code{alpha}. 
#' @param remove_last Should the last level be removed? Should be set equal to 
#' \code{TRUE} for CDF computations and \code{FALSE} for PMF computations.
#' @return List by treatment of simultaneous confidence intervals
evaluate_marg_dist_simul_ci <- function(marg_dist_est, marg_dist_eif, alpha,
                                       remove_last = FALSE){
	simul_ci <- mapply(pt_est = marg_dist_est, trt_spec_marg_dist_eif = marg_dist_eif, 
	                    FUN = compute_trt_spec_marg_dist_simul_ci,
	                    SIMPLIFY = FALSE, MoreArgs = list(remove_last = remove_last,
	                                                      alpha = alpha))
    return(simul_ci)
}

#' Compute simultaneous confidence interval for treatment-specific marginal distribution 
#' @param pt_est The point estimate of the treatment-specific marginal CDF/PMF
#' @param trt_spec_marg_dist_eif The EIF estimates for the treatment-specific marginal 
#' CDF/PMF estimates
#' @param alpha Confidence intervals have nominal level 1-\code{alpha}. 
#' @param remove_last Should the last level be removed? Should be set equal to 
#' \code{TRUE} for CDF computations and \code{FALSE} for PMF computations.
#' @importFrom stats cov quantile
#' @return Confidence interval 
compute_trt_spec_marg_dist_simul_ci <- function(pt_est, trt_spec_marg_dist_eif,
                                               remove_last = TRUE, alpha){
	# remove largest value
	if(remove_last){ # for CDF since last pt_est is always 1
		pt_est <- pt_est[-length(pt_est)] 
	}
	K <- length(pt_est)
	gradient <- diag(1 / (pt_est - pt_est^2))
	cor_mat <- stats::cor(trt_spec_marg_dist_eif %*% gradient)
	# put on logistic scale
	cov_est_logistic <- stats::cov(trt_spec_marg_dist_eif %*% gradient) / length(trt_spec_marg_dist_eif[,1])
	# Sigma <- n * cov_est_logistic
	normal_samples <- MASS::mvrnorm(n = 1e5, mu = rep(0, K),
	                                Sigma = cor_mat)
	max_samples <- apply(abs(normal_samples), 1, max)
	q_1alpha <- stats::quantile(max_samples, p = 1 - alpha)
	return(stats::plogis(stats::qlogis(pt_est) + t(c(-q_1alpha, q_1alpha) %o% sqrt(diag(cov_est_logistic)))))
}

#' Evaluate pointwise confidence interval for marginal CDF. 
#' @param marg_cdf_est The point estimate of the marginal CDF distribution
#' @param marg_cdf_eif The EIF estimates for the marginal CDF estimates
#' @param alpha Confidence intervals have nominal level 1-\code{alpha}. 
#' @return List by treatment of simultaneous confidence intervals
evaluate_marg_cdf_ptwise_ci <- function(marg_cdf_est, marg_cdf_eif, alpha){
	marg_cdf_cov <- lapply(marg_cdf_eif, function(x){
		cov(x) / length(x[,1])
	})
	# do on logistic scale
	ptwise_ci <- mapply(pt_est = marg_cdf_est, cov_est = marg_cdf_cov, 
	       				FUN = compute_trt_spec_marg_dist_ptwise_ci, 
	       				MoreArgs = list(alpha = alpha, cdf = TRUE), SIMPLIFY = FALSE)
    return(ptwise_ci)
}

#' Evaluate pointwise confidence interval for marginal PMF. 
#' @param marg_pmf_est The point estimate of the marginal PMF distribution
#' @param marg_pmf_eif The EIF estimates for the marginal PMF estimates
#' @param alpha Confidence intervals have nominal level 1-\code{alpha}. 
#' @return List by treatment of simultaneous confidence intervals
evaluate_marg_pmf_ptwise_ci <- function(marg_pmf_est, marg_pmf_eif, alpha){
	marg_pmf_cov <- lapply(marg_pmf_eif, function(x){
		cov(x) / length(x[,1])
	})
	# do on logistic scale
	ptwise_ci <- mapply(pt_est = marg_pmf_est, cov_est = marg_pmf_cov, 
	       				FUN = compute_trt_spec_marg_dist_ptwise_ci, 
	       				MoreArgs = list(alpha = alpha, cdf = FALSE), SIMPLIFY = FALSE)
    return(ptwise_ci)
}

#' Compute simultaneous confidence interval for treatment-specific marginal distribution 
#' @param pt_est The point estimate of the treatment-specific marginal CDF/PMF
#' @param cov_est Covariance matrix estimates.
#' @param alpha Confidence intervals have nominal level 1-\code{alpha}. 
#' @param cdf Is this for CDF or PMF?
#' @return Confidence interval 

compute_trt_spec_marg_dist_ptwise_ci <- function(pt_est, cov_est, alpha, cdf = TRUE){
	 K <- length(pt_est)
	 # remove largest value
	 if(cdf){
	 	pt_est <- pt_est[-K]
	 	cov_est <- cov_est[1:(K-1), 1:(K-1)]
	 }else{
	 	pt_est <- pt_est
	 	cov_est <- cov_est
	 }
	 # put on logistic scale
	 gradient <- diag(1 / (pt_est - pt_est^2))
	 cov_est_logistic <- t(gradient) %*% cov_est %*% gradient
	 return(stats::plogis(stats::qlogis(pt_est) + t(qnorm(c(alpha/2, 1 - alpha/2)) %o% sqrt(diag(cov_est_logistic)))))
}

#' Get eif estimates for treatment-specific PMF
#' 
#' @param pmf_est Estimated conditional PMF for \code{trt_level}. 
#' @param treat_prob_est Estimated propensity for \code{trt_level}.
#' @param out A \code{numeric} vector containing the outcomes. Missing outcomes are 
#' allowed. 
#' @param treat A \code{numeric} vector containing treatment status. Should only assume 
#' a value 0 or 1. 
#' @param out_levels A \code{numeric} vector containing all ordered levels of the 
#' outcome. 
#' @return a list of eif estimates
evaluate_marg_pmf_eif <- function(pmf_est, treat_prob_est, treat, out, out_levels){
	eif_matrix_list <- mapply(trt_spec_pmf_est = pmf_est, 
	       trt_spec_prob_est = treat_prob_est, trt_level = list(1,0), 
	       FUN = evaluate_trt_spec_pmf_eif,
	       MoreArgs = list(treat = treat, out = out, out_levels = out_levels),
	       SIMPLIFY = FALSE)

	return(eif_matrix_list)
}

#' Get eif estimates for treatment-specific CDF
#' 
#' @param cdf_est Estimated conditional CDF for \code{trt_level}. 
#' @param treat_prob_est Estimated propensity for \code{trt_level}.
#' @param out A \code{numeric} vector containing the outcomes. Missing outcomes are 
#' allowed. 
#' @param treat A \code{numeric} vector containing treatment status. Should only assume 
#' a value 0 or 1. 
#' @param out_levels A \code{numeric} vector containing all ordered levels of the 
#' outcome. 
#' @return a list of eif estimates

evaluate_marg_cdf_eif <- function(cdf_est, treat_prob_est, treat, out, out_levels){
	eif_matrix_list <- mapply(trt_spec_cdf_est = cdf_est, 
	       trt_spec_prob_est = treat_prob_est, trt_level = list(1,0), 
	       FUN = evaluate_trt_spec_theta_eif,
	       MoreArgs = list(treat = treat, out = out, out_levels = out_levels),
	       SIMPLIFY = FALSE)

	return(eif_matrix_list)
}

#' Marginalize over empirical distribution to obtain marginal
#' treatment-specific PMF estimate.
#' 
#' @param pmf_est Estimates of treatment-specific conditional PMF. 

marginalize_pmf <- function(pmf_est){
	lapply(pmf_est, colMeans)
}
