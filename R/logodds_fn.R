#' Compute confidence interval/s for the log-odds parameters
#' 
#' @param out A \code{numeric} vector containing the outcomes. Missing outcomes are 
#' allowed. 
#' @param treat A \code{numeric} vector containing treatment status. Should only assume 
#' a value 0 or 1. 
#' @param covar A \code{data.frame} containing the covariates to include in the working
#' proportional odds model. 
#' @param logodds_est The point estimates for log-odds.
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
#' @param cdf_est A list of treatment-specific CDF estimates.
#' @param ... Other options (not currently used). 
#' @return List with \code{wald} and \code{bca}-estimated confidence intervals 
#' for the weighted mean parameters. 
#' 
estimate_ci_logodds <- function(logodds_est, cdf_est, out_form, covar, 
                                treat_prob_est, treat, treat_form, out, ci, 
                                alpha = 0.05, nboot, out_levels, out_model, ...){
	# get ci
	if("wald" %in% ci){
		theta_cov <- evaluate_theta_cov(cdf_est = cdf_est, 
		                                treat_prob_est = treat_prob_est, 
		                                treat = treat, 
		                                out = out, out_levels = out_levels)
		beta_cov_est <- evaluate_beta_cov(cdf_est = cdf_est, 
		                                  theta_cov = theta_cov)
		# treatment 1 ci
		wald_ci_1 <- logodds_est[1] + qnorm(c(alpha/2, 1 - alpha/2)) * sqrt(beta_cov_est[1])
		# treatment 0 ci
		wald_ci_0 <- logodds_est[2] + qnorm(c(alpha/2, 1 - alpha/2)) * sqrt(beta_cov_est[2])
		# difference
		g <- matrix(c(1,-1), nrow = 2)
		wald_ci_diff <- logodds_est[3] + qnorm(c(alpha/2, 1 - alpha/2)) * sqrt(beta_cov_est[[3]])
		# format
		wald_ci <- rbind(wald_ci_1, wald_ci_0, wald_ci_diff)
	}else{
		wald_ci <- NULL
	}

	if("bca" %in% ci){
		bca_ci <- bca_logodds(treat = treat, covar = covar, 
		                      out = out, nboot = nboot, 
                      		  treat_form = treat_form, 
                      		  out_levels = out_levels, 
                      		  out_form = out_form,
                      		  logodds_est = logodds_est, 
                      		  alpha = alpha,
                      		  out_model = out_model)
	}else{
		bca_ci <- NULL
	}
	return(list(wald = wald_ci, bca = bca_ci))
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
#' @param logodds_est The estimated log-odds. 
#' @param alpha Level of confidence interval.
#' @return matrix with treatment-specific log-odds CIs and CI for difference.
bca_logodds <- function(treat, covar, out, nboot, 
                      treat_form, out_levels, out_form, out_model,
                      logodds_est, alpha = 0.05){
	boot_samples <- replicate(nboot, 
	                          one_boot_logodds(treat = treat, 
                                               covar = covar, 
	                                           out = out, 
	                                           treat_form = treat_form, 
	                                           out_levels = out_levels, 
	                                           out_form = out_form,
	                                           out_model = out_model))
	boot_trt1 <- boot_samples[1,]
	boot_trt0 <- boot_samples[2,]
	boot_diff <- boot_samples[3,]

	# remove Inf
	boot_trt1 <- boot_trt1[boot_trt1 != Inf & boot_trt1 != -Inf]
	boot_trt0 <- boot_trt0[boot_trt0 != Inf & boot_trt0 != -Inf]
	boot_diff <- boot_diff[boot_diff != Inf & boot_diff != -Inf]

	jack_samples <- jack_logodds(treat = treat,
	                           covar = covar,
	                           out = out, 
	                           treat_form = treat_form, 
	                           out_levels = out_levels, 
	                           out_form = out_form,
	                           out_model = out_model)

	jack_trt1 <- jack_samples[1,]
	jack_trt0 <- jack_samples[2,]
	jack_diff <- jack_samples[3,]
	jack_trt1 <- jack_trt1[jack_trt1 != Inf & jack_trt1 != -Inf]
	jack_trt0 <- jack_trt0[jack_trt0 != Inf & jack_trt0 != -Inf]
	jack_diff <- jack_diff[jack_diff != Inf & jack_diff != -Inf]
	
	# CI for 1 
	bca_ci_trt1 <- bca_interval(pt_est = logodds_est[1],
	                            boot_samples = boot_trt1,
	                            jack_samples = jack_trt1,
	                            alpha = alpha)

	# CI for 0 
	bca_ci_trt0 <- bca_interval(pt_est = logodds_est[2],
                            boot_samples = boot_trt0,
                            jack_samples = jack_trt0,
                            alpha = alpha)

	# CI for diff
	bca_ci_diff <- bca_interval(pt_est = logodds_est[3],
                        boot_samples = boot_diff,
                        jack_samples = jack_diff,
                        alpha = alpha)
	return(rbind(bca_ci_trt1, bca_ci_trt0, bca_ci_diff))
}

#' Compute jackknife log-odds estimates.
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
#' @return Jackknife estimated log-odds
jack_logodds <- function(treat, covar, out, treat_form, out_model, out_levels, out_form){
  	logodds_jack_est <- sapply(seq_along(out), function(i){
		logodds_minusi <- get_one_logodds(treat = treat[-i],
		                                covar = covar[-i, , drop = FALSE],
		                                out = out[-i],
		                                treat_form = treat_form,
		                                out_levels = out_levels,
		                                out_form = out_form,
		                                out_model = out_model)  		
		return(logodds_minusi)
  	})
  	return(logodds_jack_est)
}

#' Get one bootstrap computation of the log odds parameters. 
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
#' @return Estimates of log odds for a particular bootstrap sample. 
one_boot_logodds <- function(treat, covar, out, treat_form, 
                             out_levels, out_form, out_model){
	boot_idx <- sample(seq_along(out), replace = TRUE)
	logodds_boot_est <- tryCatch({get_one_logodds(treat = treat[boot_idx],
	                                covar = covar[boot_idx, , drop = FALSE],
	                                out = out[boot_idx],
	                                treat_form = treat_form,
	                                out_levels = out_levels,
	                                out_form = out_form,
	                                out_model = out_model)}, error = function(e){
		rep(NA, 3)
	})
	return(logodds_boot_est)
}

#' Compute one log odds based on a given data set. 
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
#' @return Estimated log odds for these input data. 
get_one_logodds <- function(treat, covar, treat_form, out_model,
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
	                        out_model = out_model,
	                        out_form = out_form, 
	                        treat_prob_est = treat_prob_est,
	                        return_models = FALSE)
  	pmf_est <- pmf_fit$pmf

	cdf_est <- estimate_cdf(pmf_est = pmf_est)

  	logodds_est <- estimate_logodds(cdf_est = cdf_est)

  	return(logodds_est)
}

#' implements a plug-in estimator of equation (2) in Diaz et al
#' @param cdf_est A list of treatment-specific CDF estimates
#' @return Log odds of treatment = 1, = 0, and the difference.  
estimate_logodds <- function(cdf_est){
	# get marginal CDF
	theta_1 <- colMeans(cdf_est[[1]])
	theta_0 <- colMeans(cdf_est[[2]])

	# projection
	K <- length(theta_1)
	beta_1 <- mean(qlogis(theta_1[1:(K-1)]))
	beta_0 <- mean(qlogis(theta_0[1:(K-1)]))
	beta_est <- beta_1 - beta_0
	return(c(beta_1, beta_0, beta_est))
}

#' Get the covariance matrix for beta
#' @param cdf_est Estimated CDFs
#' @param theta_cov Covariance matrix for CDF estimates
#' @return Estimated covariance matrix for log-odds ratio parameters 
evaluate_beta_cov <- function(cdf_est, theta_cov){
	# get marginal CDF
	theta_1 <- colMeans(cdf_est[[1]])
	theta_0 <- colMeans(cdf_est[[2]])
	K <- length(theta_1)
	avg_vec <- rep(1 / (K-1), K-1)
	zero_vec <- rep(0, K-1)
	# gradient
	grad_1 <- avg_vec * 
			sapply(theta_1[1:(K-1)], function(x){
				1 / (x - x^2) # deriv of log(x / (1 - x))
			})
	grad_0 <- avg_vec * 
				sapply(theta_0[1:(K-1)], function(x){
					1 / (x - x^2) # deriv of log(x / (1 - x))
				})
	grad_diff <- c(avg_vec, -avg_vec) * 
				sapply(c(theta_1[1:(K-1)], theta_0[1:(K-1)]), function(x){
					1 / (x - x^2) # deriv of log(x / (1 - x))
				})
	theta_1_avg_cov_est <- t(c(grad_1, zero_vec)) %*% theta_cov %*% c(grad_1, zero_vec)
	theta_0_avg_cov_est <- t(c(zero_vec, grad_0)) %*% theta_cov %*% c(zero_vec, grad_0)
	beta_cov_est <- t(grad_diff) %*% theta_cov %*% grad_diff
	return(c(theta_1_avg_cov_est, theta_0_avg_cov_est, beta_cov_est))
}

#' get a covariance matrix for the estimated CDF
#' @param cdf_est The estimates of the treatment-specific CDFs
#' @param treat_prob_est List of estimated probability of treatments, output from call
#' to \code{estimate_treat_prob}.
#' @param out A \code{numeric} vector containing the outcomes. Missing outcomes are 
#' allowed. 
#' @param treat A \code{numeric} vector containing treatment status. Should only assume 
#' a value 0 or 1. 
#' @param out_levels A \code{numeric} vector containing all ordered levels of the 
#' outcome. 
#' @return Estimated covariance matrix for CDF estimates
evaluate_theta_cov <- function(cdf_est, treat_prob_est, treat, out, out_levels){
	eif_matrix_list <- mapply(trt_spec_cdf_est = cdf_est, 
	       trt_spec_prob_est = treat_prob_est, trt_level = list(1,0), 
	       FUN = evaluate_trt_spec_theta_eif,
	       MoreArgs = list(treat = treat, out = out, out_levels = out_levels),
	       SIMPLIFY = FALSE)
	eif_matrix <- Reduce(cbind, eif_matrix_list)
	cov_matrix <- cov(eif_matrix) / length(out)
	return(cov_matrix)
}

#' get a matrix of eif estimates for the treatment-specific CDF estimates
#' @param trt_spec_cdf_est Estimated conditional CDF for \code{trt_level}. 
#' @param trt_spec_prob_est Estimated propensity for \code{trt_level}.
#' @param out A \code{numeric} vector containing the outcomes. Missing outcomes are 
#' allowed. 
#' @param treat A \code{numeric} vector containing treatment status. Missing
#' values are not allowed unless the corresponding entry in \code{out} is also missing. 
#' Only values of 0 or 1 are treated as actual treatment levels. Any other value is assumed 
#' to encode a value for which the outcome is missing and the corresponding outcome value is 
#' ignored. 
#' @param trt_level Treatment level 
#' @param out_levels A \code{numeric} vector containing all ordered levels of the 
#' outcome. 
#' @return matrix of EIF estimates for CDF. 
evaluate_trt_spec_theta_eif <- function(trt_spec_cdf_est, 
                               trt_spec_prob_est, 
                               trt_level,
                               treat, out, out_levels){
	K <- ncol(trt_spec_cdf_est)
	eif_matrix <- matrix(NA, nrow = length(out), ncol = ncol(trt_spec_cdf_est) - 1)
	for(k in 1:(K-1)){
		eif_matrix[,k] <- eif_theta_k(k = out_levels[k], out = out, treat = treat, 
		                              trt_level = trt_level, 
		                              trt_spec_prob_est = trt_spec_prob_est,
		                              trt_k_spec_cdf_est = trt_spec_cdf_est[,k])
	}
	return(eif_matrix)
}

#' Get a matrix of eif estimates for treatment-specific PMF
#' 
#' @param trt_spec_pmf_est Estimated conditional PMF for \code{trt_level}. 
#' @param trt_spec_prob_est Estimated propensity for \code{trt_level}.
#' @param out A \code{numeric} vector containing the outcomes. Missing outcomes are 
#' allowed. 
#' @param treat A \code{numeric} vector containing treatment status. Should only assume 
#' a value 0 or 1. 
#' @param trt_level Treatment level 
#' @param out_levels A \code{numeric} vector containing all ordered levels of the 
#' outcome. 
#' @return a matrix of EIF estimates
evaluate_trt_spec_pmf_eif <- function(trt_spec_pmf_est, 
                               trt_spec_prob_est, 
                               trt_level,
                               treat, out, out_levels){
	K <- ncol(trt_spec_pmf_est)
	eif_matrix <- matrix(NA, nrow = length(out), ncol = ncol(trt_spec_pmf_est))
	for(k in seq_len(K)){
		eif_matrix[,k] <- eif_pmf_k(k = out_levels[k], out = out, treat = treat, 
		                              trt_level = trt_level, 
		                              trt_spec_prob_est = trt_spec_prob_est,
		                              trt_k_spec_pmf_est = trt_spec_pmf_est[,k])
	}
	return(eif_matrix)
}

#' Get EIF estimates for treatment-specific PMF at a particular
#' level of the outcome
#' 
#' @param k The level of the outcome.
#' @param trt_k_spec_pmf_est Estimated conditional PMF for \code{trt_level} at \code{k}. 
#' @param trt_spec_prob_est Estimated propensity for \code{trt_level}.
#' @param out A \code{numeric} vector containing the outcomes. Missing outcomes are 
#' allowed. 
#' @param treat A \code{numeric} vector containing treatment status. Should only assume 
#' a value 0 or 1. 
#' @param trt_level Treatment level 

eif_pmf_k <- function(k, out, treat, trt_level, trt_spec_prob_est,
                        trt_k_spec_pmf_est){
	out[is.na(out)] <- -99999
	eif <- as.numeric(treat == trt_level) / trt_spec_prob_est * 
		(as.numeric(out == k) - trt_k_spec_pmf_est) + 
			trt_k_spec_pmf_est - mean(trt_k_spec_pmf_est)
	return(eif)
}


#' Get EIF estimates for treatment-specific CDF at a particular
#' level of the outcome
#' 
#' @param k The level of the outcome.
#' @param trt_k_spec_cdf_est Estimated conditional CDF for \code{trt_level} at \code{k}. 
#' @param trt_spec_prob_est Estimated propensity for \code{trt_level}.
#' @param out A \code{numeric} vector containing the outcomes. Missing outcomes are 
#' allowed. 
#' @param treat A \code{numeric} vector containing treatment status. Should only assume 
#' a value 0 or 1. 
#' @param trt_level Treatment level 

eif_theta_k <- function(k, out, treat, trt_level, trt_spec_prob_est,
                        trt_k_spec_cdf_est){
	out[is.na(out)] <- -99999
	eif <- as.numeric(treat == trt_level) / trt_spec_prob_est * 
		(as.numeric(out <= k) - trt_k_spec_cdf_est) + 
			trt_k_spec_cdf_est - mean(trt_k_spec_cdf_est)
	return(eif)
}

