#' Compute the estimate of the weighted mean parameter based on 
#' estimated PMF in each treatment arm.  
#' 
#' @param pmf_est List of treatment-specific PMF estimates.
#' @param out A \code{numeric} vector containing the outcomes. Missing outcomes are 
#' allowed. 
#' @param treat A \code{numeric} vector containing treatment status. Missing
#' values are not allowed unless the corresponding entry in \code{out} is also missing. 
#' Only values of 0 or 1 are treated as actual treatment levels. Any other value is assumed 
#' to encode a value for which the outcome is missing and the corresponding outcome value is 
#' ignored. 
#' @param out_levels A \code{numeric} vector containing all ordered levels of the 
#' outcome. 
#' @param out_weights A vector of \code{numeric} weights with length equal to the length 
#' of \code{out_levels}. 
#' @param treat_prob_est Estimated probability of treatments, output from call
#' to \code{estimate_treat_prob}.
#' @param return_cov If \code{TRUE} the estimated covariance matrix is returned.
#' @return List with estimates of treatment-specific means and difference in means. 
#' If \code{return_cov = TRUE}, also includes covariance matrix estimates. 
estimate_wmean <- function(pmf_est, treat, out, out_levels, out_weights,
                           treat_prob_est, return_cov = TRUE){
	# align ordering 
	ordered_out_levels <- out_levels[order(out_levels)]
	ordered_out_weights <- out_weights[order(out_weights)]

	# map into treatment-specific conditional mean estimate
	cond_mean_est <- lapply(pmf_est, estimate_cond_mean, 
	                        ordered_out_levels = ordered_out_levels,
	                        ordered_out_weights = ordered_out_weights)

	# influence function
	wmean_eif_est <- mapply(trt_spec_cond_mean_est = cond_mean_est,
	                        trt_spec_prob_est = treat_prob_est,
	                        trt_level = 1:0, 
	                        FUN = estimate_eif_wmean,
	                        MoreArgs = list(out = out, treat = treat),
	                        SIMPLIFY = FALSE)
	# point estimate
	wmean_est <- sapply(cond_mean_est, mean)
	wmean_est <- c(wmean_est, wmean_est[1] - wmean_est[2])
	rslt <- list(est = wmean_est, cov = NULL)

	if(return_cov){
		eif_mat <- Reduce(cbind, wmean_eif_est)
		colnames(eif_mat) <- NULL
		rslt$cov <- cov(eif_mat) / length(out)
	}
	return(rslt)
}

#' Map an estimate of treatment-specific PMF into an estimate of 
#' treatment specific conditional mean for each observation.
#' @param trt_spec_pmf_est The treatment-specific PMF estimates 
#' @param ordered_out_levels Self explanatory
#' @param ordered_out_weights Self explanatory
#' @return Vector of estimated conditional means
estimate_cond_mean <- function(trt_spec_pmf_est, ordered_out_levels, ordered_out_weights){
	rslt <- apply(trt_spec_pmf_est, 1, function(pmf){
		sum(ordered_out_levels * ordered_out_weights * pmf) / sum(ordered_out_weights * pmf)
	})
	return(rslt)
}

#' Obtain an estimate of the efficient influence function for the 
#' treatment-specific weighted mean parameter
#' @param trt_spec_cond_mean_est Conditional mean for \code{trt_level}
#' @param trt_spec_prob_est Propensity for \code{trt_level}
#' @param trt_level Treatment level 
#' @param out A \code{numeric} vector containing the outcomes. Missing outcomes are 
#' allowed. 
#' @param treat A \code{numeric} vector containing treatment status. Missing
#' values are not allowed unless the corresponding entry in \code{out} is also missing. 
#' Only values of 0 or 1 are treated as actual treatment levels. Any other value is assumed 
#' to encode a value for which the outcome is missing and the corresponding outcome value is 
#' ignored. 
estimate_eif_wmean <- function(trt_spec_cond_mean_est,
                              trt_spec_prob_est, 
                              trt_level,
                              out, treat){
	out[is.na(out)] <- -99999
	eif <- as.numeric(treat == trt_level) / trt_spec_prob_est * (out - trt_spec_cond_mean_est) + 
			trt_spec_cond_mean_est - mean(trt_spec_cond_mean_est)
	return(eif)
}

#' Compute confidence interval/s for the weight mean parameters
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
#' @param wmean_est The point estimates for weighted means
#' @param alpha Confidence intervals have nominal level 1-\code{alpha}. 
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
#' @param ci A vector of \code{characters} indicating which confidence intervals
#' should be computed (\code{"bca"} and/or \code{"wald"}) 
#' @param nboot Number of bootstrap replicates used to compute bootstrap confidence
#' intervals. 
#' @return List with \code{wald} and \code{bca}-estimated confidence intervals 
#' for the weighted mean parameters. 
estimate_ci_wmean <- function(
  out,
  treat,
  covar,
  wmean_est,
  alpha = 0.05, 
  out_levels = order(unique(out)),
  out_form = NULL,
  out_weights = rep(1, length(out_levels)),
  out_model,
  treat_form = "1",
  ci = c("bca", "wald"), 
  nboot = 1e4
  ){
	if("wald" %in% ci){
		wald_ci <- wald_ci_wmean(wmean_est, alpha = alpha)
	}else{
		wald_ci <- NULL
	}

	if("bca" %in% ci){
		bca_ci <- bca_wmean(treat = treat, covar = covar,
		                    out = out, nboot = nboot, 
		                    treat_form = treat_form, out_levels = out_levels,
		                    out_form = out_form, wmean_est = wmean_est, 
		                    alpha = alpha, out_weights = out_weights,
		                    out_model = out_model)
	}else{
		bca_ci <- NULL
	}
	return(list(wald = wald_ci, bca = bca_ci))
}

#' Compute a Wald confidence interval for the weighted mean
#' @param wmean_est The estimated weighted means + estimated covariance matrix. 
#' @param alpha Level of confidence interval.
#' @return matrix with treatment-specific weighted mean CIs and CI for difference.
wald_ci_wmean <- function(wmean_est, alpha){
	# treatment 1 ci
	wald_ci_1 <- wmean_est$est[1] + qnorm(c(alpha/2, 1 - alpha/2)) * sqrt(wmean_est$cov[1,1])
	# treatment 0 ci
	wald_ci_0 <- wmean_est$est[2] + qnorm(c(alpha/2, 1 - alpha/2)) * sqrt(wmean_est$cov[2,2])
	# difference
	g <- matrix(c(1,-1), nrow = 2)
	wald_ci_diff <- wmean_est$est[1] - wmean_est$est[2] + qnorm(c(alpha/2, 1 - alpha/2)) * c(sqrt(t(g)%*%wmean_est$cov%*%g))
	# format
	return(rbind(wald_ci_1,wald_ci_0,wald_ci_diff))
}

#' Compute a BCa bootstrap confidence interval for the weighted mean. The code is 
#' based on the slides found here: http://users.stat.umn.edu/~helwig/notes/bootci-Notes.pdf
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
#' @param out_weights A vector of \code{numeric} weights with length equal to the length 
#' of \code{out_levels}. 
#' @param wmean_est The estimated weighted means + estimated covariance matrix. 
#' @param alpha Level of confidence interval.
#' @return matrix with treatment-specific weighted mean CIs and CI for difference.
bca_wmean <- function(treat, covar, out, nboot, 
                      treat_form, out_levels, out_form, out_weights,
                      out_model, 
                      wmean_est, alpha = 0.05){
	boot_samples <- replicate(nboot, 
	                          one_boot_wmean(treat = treat, 
                                             covar = covar, 
	                                         out = out, 
	                                         treat_form = treat_form, 
	                                         out_levels = out_levels, 
	                                         out_form = out_form,
	                                         out_weights = out_weights,
	                                         out_model = out_model))
	boot_trt1 <- boot_samples[1,]
	boot_trt0 <- boot_samples[2,]
	boot_diff <- boot_samples[1,] - boot_samples[2,]

	jack_samples <- jack_wmean(treat = treat,
	                           covar = covar,
	                           out = out, 
	                           treat_form = treat_form, 
	                           out_levels = out_levels, 
	                           out_form = out_form,
	                           out_weights = out_weights,
	                           out_model = out_model)

	jack_trt1 <- jack_samples[1,]
	jack_trt0 <- jack_samples[2,]
	jack_diff <- jack_samples[1,] - jack_samples[2,]

	# CI for 1 
	bca_ci_trt1 <- bca_interval(pt_est = wmean_est$est[1],
	                            boot_samples = boot_trt1,
	                            jack_samples = jack_trt1,
	                            alpha = alpha)

	# CI for 0 
	bca_ci_trt0 <- bca_interval(pt_est = wmean_est$est[2],
                            boot_samples = boot_trt0,
                            jack_samples = jack_trt0,
                            alpha = alpha)

	# CI for diff
	bca_ci_diff <- bca_interval(pt_est = wmean_est$est[1] - wmean_est$est[2],
                        boot_samples = boot_diff,
                        jack_samples = jack_diff,
                        alpha = alpha)
	return(rbind(bca_ci_trt1, bca_ci_trt0, bca_ci_diff))
}

#' Compute jackknife weighted mean estimates.
#' @param out A \code{numeric} vector containing the outcomes. Missing outcomes are 
#' allowed. 
#' @param treat A \code{numeric} vector containing treatment status. Missing
#' values are not allowed unless the corresponding entry in \code{out} is also missing. 
#' Only values of 0 or 1 are treated as actual treatment levels. Any other value is assumed 
#' to encode a value for which the outcome is missing and the corresponding outcome value is 
#' ignored. 
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
#' @param out_weights A vector of \code{numeric} weights with length equal to the length 
#' of \code{out_levels}. 
#' @return Jackknife-estimated weighted mean

jack_wmean <- function(treat, covar, out, treat_form, out_levels, 
                       out_form, out_weights, out_model){
  	wmean_jack_est <- sapply(seq_along(out), function(i){
		wmean_minusi <- get_one_wmean(treat = treat[-i],
		                                covar = covar[-i, , drop = FALSE],
		                                out = out[-i],
		                                treat_form = treat_form,
		                                out_levels = out_levels,
		                                out_form = out_form,
		                                out_weights = out_weights,
		                                out_model = out_model)  		
		return(wmean_minusi)
  	})
  	return(wmean_jack_est)
}

#' Get one bootstrap computation of the weighted mean parameters. 
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
#' @param treat_form The right-hand side of a regression formula for the working model of
#' treatment probability as a function of covariates
#' @param out_levels A \code{numeric} vector containing all ordered levels of the 
#' outcome. 
#' @param out_form The right-hand side of a regression formula for the working proportional 
#' odds model. NOTE: THIS FORMULA MUST NOT SUPPRESS THE INTERCEPT. 
#' @param out_model Which R function should be used to fit the proportional odds 
#' model. Options are \code{"polr"} (from the \code{MASS} package), 
#' "vglm" (from the \code{VGAM} package), or \code{"clm"} (from the \code{ordinal} package).
#' @param out_weights A vector of \code{numeric} weights with length equal to the length 
#' of \code{out_levels}. 
#' @return Estimates of weighted mean for a particular bootstrap sample. 
one_boot_wmean <- function(treat, covar, out, treat_form, out_levels, 
                           out_form, out_weights, out_model){
	boot_idx <- sample(seq_along(out), replace = TRUE)
	wmean_boot_est <- tryCatch({get_one_wmean(treat = treat[boot_idx],
	                                covar = covar[boot_idx, , drop = FALSE],
	                                out = out[boot_idx],
	                                treat_form = treat_form,
	                                out_levels = out_levels,
	                                out_form = out_form,
	                                out_weights = out_weights,
	                                out_model = out_model)}, error = function(e){
		rep(NA, 3)
	})
	return(wmean_boot_est)
}

#' Compute one weighted mean based on a given data set. 
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
#' @param treat_form The right-hand side of a regression formula for the working model of
#' treatment probability as a function of covariates
#' @param out_levels A \code{numeric} vector containing all ordered levels of the 
#' outcome. 
#' @param out_form The right-hand side of a regression formula for the working proportional 
#' odds model. NOTE: THIS FORMULA MUST NOT SUPPRESS THE INTERCEPT. 
#' @param out_model Which R function should be used to fit the proportional odds 
#' model. Options are \code{"polr"} (from the \code{MASS} package), 
#' "vglm" (from the \code{VGAM} package), or \code{"clm"} (from the \code{ordinal} package).
#' @param out_weights A vector of \code{numeric} weights with length equal to the length 
#' of \code{out_levels}. 
get_one_wmean <- function(treat, covar, treat_form,
                          out, out_levels, out_form,
                          out_model,
                          out_weights){
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
	                        out_model = out_model, return_models = FALSE)
  	pmf_est <- pmf_fit$pmf

	
  	wmean_est <- estimate_wmean(
          pmf_est = pmf_est, treat = treat, out = out, out_levels = out_levels, 
          out_weights = out_weights, treat_prob_est = treat_prob_est, 
          return_cov = FALSE
        )
  	return(wmean_est$est)
}