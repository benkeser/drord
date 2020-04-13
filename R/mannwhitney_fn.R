
estimate_ci_mannwhitney <- function(
    mannwhitney_est, cdf_est, pmf_est, treat_prob_est, treat_form, out_form,
    treat, ci, out, alpha, nboot, out_levels, covar
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
                      		  alpha = alpha)
	}else{
		bca_ci <- NULL
	}

	return(list(wald = wald_ci, bca = bca_ci))
}

bca_mannwhitney <- function(treat, covar, out, nboot, 
                      treat_form, out_levels, out_form,
                      mannwhitney_est, alpha = 0.05){
	boot_samples <- replicate(nboot, 
	                          one_boot_mannwhitney(treat = treat, 
                                               covar = covar, 
	                                           out = out, 
	                                           treat_form = treat_form, 
	                                           out_levels = out_levels, 
	                                           out_form = out_form))

	jack_samples <- jack_mannwhitney(treat = treat,
	                           covar = covar,
	                           out = out, 
	                           treat_form = treat_form, 
	                           out_levels = out_levels, 
	                           out_form = out_form)

	# CI for 1 
	bca_ci_mannwhitney <- bca_interval(pt_est = mannwhitney_est,
	                            boot_samples = boot_samples,
	                            jack_samples = jack_samples,
	                            alpha = alpha)

	return(rbind(bca_ci_mannwhitney))
}

jack_mannwhitney <- function(treat, covar, out, treat_form, out_levels, out_form){
  	mannwhitney_jack_est <- sapply(seq_along(out), function(i){
		mannwhitney_minusi <- get_one_mannwhitney(treat = treat[-i],
		                                covar = covar[-i, , drop = FALSE],
		                                out = out[-i],
		                                treat_form = treat_form,
		                                out_levels = out_levels,
		                                out_form = out_form)  		
		return(mannwhitney_minusi)
  	})
  	return(mannwhitney_jack_est)
}

one_boot_mannwhitney <- function(treat, covar, out, treat_form, out_levels, out_form){
	boot_idx <- sample(seq_along(out), replace = TRUE)
	mannwhitney_boot_est <- tryCatch({get_one_mannwhitney(treat = treat[boot_idx],
	                                covar = covar[boot_idx, , drop = FALSE],
	                                out = out[boot_idx],
	                                treat_form = treat_form,
	                                out_levels = out_levels,
	                                out_form = out_form)}, error = function(e){
		NA
	})
	return(mannwhitney_boot_est)
}

get_one_mannwhitney <- function(treat, covar, treat_form,
                            	out, out_levels, out_form){
	# obtain estimate of treatment probabilities
	treat_prob_est <- estimate_treat_prob(treat = treat,
	                                      covar = covar,
	                                      treat_form = treat_form)

	# obtain estimate of conditional PMF for each treatment level
	pmf_est <- estimate_pmf(out = out, treat = treat, 
	                        covar = covar, out_levels = out_levels,
	                        out_form = out_form, treat_prob_est = treat_prob_est)

	cdf_est <- estimate_cdf(pmf_est = pmf_est)

  	mannwhitney_est <- estimate_mannwhitney(cdf_est = cdf_est, pmf_est = pmf_est)

  	return(mannwhitney_est)
}

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

estimate_mannwhitney <- function(cdf_est, pmf_est){
	K <- ncol(cdf_est[[1]])

	# get marginal CDF
	F_1 <- colMeans(cdf_est[[1]])
	F_0 <- colMeans(cdf_est[[2]])
	
	# get marginal PDF
	f_1 <- colMeans(pmf_est[[1]])
	f_0 <- colMeans(pmf_est[[1]])

	# F(k-1 | A = 0), k = 0, ..., K
	F_0_kminus1 <- c(0, F_0[-K])

	# estimate
	mannwhitney_est <- sum((F_0_kminus1 + 1/2 * f_0) * f_1)

	return(mannwhitney_est)
}