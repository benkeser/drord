#' Compute confidence interval/s for the Mann-Whitney parameter
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
#' @param mannwhitney_est The point estimates for log-odds.
#' @param pmf_est The estimated conditional PMF.
#' @param cdf_est The estimated conditional CDF.
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
#' @return List with \code{wald} and \code{bca}-estimated confidence intervals 
#' for the Mann-Whitney parameter. 
#' 
estimate_ci_mannwhitney <- function(
    mannwhitney_est, cdf_est, pmf_est, treat_prob_est, treat_form, out_form,
    treat, ci, out, alpha, nboot, out_levels, covar, out_model
){
	n <- length(out)
	if("wald" %in% ci){
		# get big influence function matrix
		eif_F_0 <- evaluate_trt_spec_theta_eif(
			trt_spec_cdf_est = cdf_est[[2]], 
	       	trt_spec_prob_est = treat_prob_est[[2]], 
	       	trt_level = 0, treat = treat, out = out,
	       	out_levels = out_levels
        )

		eif_f_01_list <- mapply(trt_spec_pmf_est = pmf_est, 
	       trt_spec_prob_est = treat_prob_est, trt_level = list(1,0), 
	       FUN = evaluate_trt_spec_pmf_eif,
	       MoreArgs = list(treat = treat, out = out, out_levels = out_levels),
	       SIMPLIFY = FALSE)
		eif_f_01 <- Reduce(cbind, eif_f_01_list)
		eif_matrix <- cbind(eif_F_0, eif_f_01)		
		cov_matrix <- cov(eif_matrix)
		gradient <- evaluate_mannwhitney_gradient(cdf_est = cdf_est, pmf_est = pmf_est)
		se_est_mannwhitney_est <- sqrt(t(gradient) %*% cov_matrix %*% gradient / n)
		wald_ci <- mannwhitney_est + qnorm(c(alpha/2, 1 - alpha/2)) * c(se_est_mannwhitney_est)
	}else{
		wald_ci <- NULL
	}

	if("bca" %in% ci){
		bca_ci <- bca_mannwhitney(treat = treat, covar = covar, 
		                      out = out, nboot = nboot, 
                      		  treat_form = treat_form, 
                      		  out_levels = out_levels, 
                      		  out_form = out_form,
                      		  mannwhitney_est = mannwhitney_est,
                      		  alpha = alpha,
                      		  out_model = out_model)
	}else{
		bca_ci <- NULL
	}

	return(list(wald = wald_ci, bca = bca_ci))
}

#' Compute a BCa bootstrap confidence interval for the Mann-Whitney parameter. The code is 
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
#' @param mannwhitney_est The point estimate of the Mann-Whitney parameter. 
#' @param alpha Level of confidence interval.
#' @return Confidence interval for the Mann-Whitney parameter
bca_mannwhitney <- function(treat, covar, out, nboot, 
                      treat_form, out_levels, out_form,
                      mannwhitney_est, 
                      out_model, alpha = 0.05){
	boot_samples <- replicate(nboot, 
	                          one_boot_mannwhitney(treat = treat, 
                                               covar = covar, 
	                                           out = out, 
	                                           treat_form = treat_form, 
	                                           out_levels = out_levels, 
	                                           out_form = out_form,
	                                           out_model = out_model))

	jack_samples <- jack_mannwhitney(treat = treat,
	                           covar = covar,
	                           out = out, 
	                           treat_form = treat_form, 
	                           out_levels = out_levels, 
	                           out_form = out_form,
	                           out_model = out_model)

	bca_ci_mannwhitney <- bca_interval(pt_est = mannwhitney_est,
	                            boot_samples = boot_samples,
	                            jack_samples = jack_samples,
	                            alpha = alpha)

	return(rbind(bca_ci_mannwhitney))
}


#' Compute Mann-Whitney log-odds estimates.
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
#' @return Jackknife estimate of Mann-Whitney parameter

jack_mannwhitney <- function(treat, covar, out, treat_form, out_levels, out_form,
                             out_model){
  	mannwhitney_jack_est <- sapply(seq_along(out), function(i){
		mannwhitney_minusi <- get_one_mannwhitney(treat = treat[-i],
		                                covar = covar[-i, , drop = FALSE],
		                                out = out[-i],
		                                treat_form = treat_form,
		                                out_levels = out_levels,
		                                out_form = out_form,
		                                out_model = out_model)  		
		return(mannwhitney_minusi)
  	})
  	return(mannwhitney_jack_est)
}

#' Get one bootstrap computation of the Mann-Whitney parameter. 
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
#' @return Estimates of Mann-Whitney parameter for a particular bootstrap sample. 
one_boot_mannwhitney <- function(treat, covar, out, treat_form, out_levels, out_form,
                                 out_model){
	boot_idx <- sample(seq_along(out), replace = TRUE)
	mannwhitney_boot_est <- tryCatch({get_one_mannwhitney(treat = treat[boot_idx],
	                                covar = covar[boot_idx, , drop = FALSE],
	                                out = out[boot_idx],
	                                treat_form = treat_form,
	                                out_levels = out_levels,
	                                out_model = out_model,
	                                out_form = out_form)}, error = function(e){
		NA
	})
	return(mannwhitney_boot_est)
}

#' Compute one estimate of Mann-Whitney parameter based on a given data set. 
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
#' @return Estimate of Mann-Whitney parameter for these input data. 

get_one_mannwhitney <- function(treat, covar, treat_form,
                            	out, out_levels, out_form,
                            	out_model){
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

  	mannwhitney_est <- estimate_mannwhitney(cdf_est = cdf_est, pmf_est = pmf_est)

  	return(mannwhitney_est)
}

#' Compute the estimated gradient of the Mann-Whitney parameter. Needed to derive 
#' standard error for Wald confidence intervals.
#' @param cdf_est Conditional CDF estimates
#' @param pmf_est Conditional PMF estimates
#' @return 3-length vector for delta method calculus
evaluate_mannwhitney_gradient <- function(cdf_est, pmf_est){
	K <- ncol(cdf_est[[1]])

	# get marginal CDF
	F_1 <- colMeans(cdf_est[[1]])
	F_0 <- colMeans(cdf_est[[2]])
	
	# get marginal PDF
	f_1 <- colMeans(pmf_est[[1]])
	f_0 <- colMeans(pmf_est[[1]])

	# F(k-1 | A = 0), k = 0, ..., K
	F_0_kminus1 <- c(0, F_0[-K])

	gradient <- c(
      f_1[-1], 
      F_0_kminus1 + 1/2 * f_0,
      1/2 * f_1
    )
    return(gradient)
}

#' Compute the estimate of Mann-Whitney based on conditional CDF and PMF
#' @param cdf_est Conditional CDF estimates
#' @param pmf_est Conditional PMF estimates
#' @return Mann-Whitney point estimate
estimate_mannwhitney <- function(cdf_est, pmf_est){
	K <- ncol(cdf_est[[1]])

	# get marginal CDF
	F_1 <- colMeans(cdf_est[[1]])
	F_0 <- colMeans(cdf_est[[2]])
	
	# get marginal PDF
	f_1 <- colMeans(pmf_est[[1]])
	f_0 <- colMeans(pmf_est[[2]])

	# F(k-1 | A = 0), k = 0, ..., K
	F_0_kminus1 <- c(0, F_0[-K])

	# estimate
	mannwhitney_est <- sum((F_0_kminus1 + 1/2 * f_0) * f_1)

	return(mannwhitney_est)
}