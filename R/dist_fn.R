
estimate_ci_marg_dist <- function(marg_cdf_est,
                                  marg_pmf_est,                                  
                                  cdf_est,
                                  pmf_est,
	                              treat_prob_est, 
	                              treat_form, out_form,
	                              treat, ci, out_levels, 
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
		marg_cdf_simul_wald_ci <- evaluate_marg_cdf_simul_ci(marg_cdf_est = marg_cdf_est,
	                                                    marg_cdf_eif = marg_cdf_eif,
	                                                    alpha = alpha, n = length(out), 
	                                                    remove_last = TRUE)

		marg_pmf_eif <- evaluate_marg_pmf_eif(pmf_est = pmf_est, 
		                                   treat_prob_est = treat_prob_est, 
		                                   treat = treat, out = out,
		                                   out_levels = out_levels)

		marg_pmf_ptwise_wald_ci <- evaluate_marg_pmf_ptwise_ci(marg_pmf_est = marg_pmf_est,
		                                                  marg_pmf_eif = marg_pmf_eif,
		                                                  alpha = alpha)
		# simultaneous CI
		marg_pmf_simul_wald_ci <- evaluate_marg_cdf_simul_ci(marg_cdf_est = marg_pmf_est,
	                                                    marg_cdf_eif = marg_pmf_eif,
	                                                    alpha = alpha, n = length(out),
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
		                       alpha = alpha)
		bca_ci_cdf <- bca_ci$cdf
		bca_ci_pmf <- bca_ci$pmf
	}else{
		bca_ci_cdf <- NULL
		bca_ci_pmf <- NULL
	}
	return(list(cdf = list(wald = wald_ci_cdf, bca = bca_ci_cdf),
	       		pmf = list(wald = wald_ci_pmf, bca = bca_ci_pmf)))
}


marginalize_cdf <- function(cdf_est){
	lapply(cdf_est, colMeans)
}

bca_marg_dist <- function(treat, covar, out, nboot, 
                      treat_form, out_levels, out_form,
                      marg_cdf_est, marg_pmf_est, alpha = 0.05){
	K <- length(out_levels)
	boot_samples <- replicate(nboot, 
	                          one_boot_marg_cdf(treat = treat, 
                                               covar = covar, 
	                                           out = out, 
	                                           treat_form = treat_form, 
	                                           out_levels = out_levels, 
	                                           out_form = out_form))

	jack_samples <- jack_marg_cdf(treat = treat,
	                           covar = covar,
	                           out = out, 
	                           treat_form = treat_form, 
	                           out_levels = out_levels, 
	                           out_form = out_form)

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

jack_marg_cdf <- function(treat, covar, out, treat_form, out_levels, out_form){
  	marg_cdf_jack_est <- sapply(seq_along(out), function(i){
		marg_cdf_minusi <- get_one_marg_cdf(treat = treat[-i],
		                                covar = covar[-i, , drop = FALSE],
		                                out = out[-i],
		                                treat_form = treat_form,
		                                out_levels = out_levels,
		                                out_form = out_form)  		
		return(marg_cdf_minusi)
  	})
  	return(marg_cdf_jack_est)
}

one_boot_marg_cdf <- function(treat, covar, out, treat_form, out_levels, out_form){
	boot_idx <- sample(seq_along(out), replace = TRUE)
	marg_cdf_boot_est <- get_one_marg_cdf(treat = treat[boot_idx],
	                                covar = covar[boot_idx, , drop = FALSE],
	                                out = out[boot_idx],
	                                treat_form = treat_form,
	                                out_levels = out_levels,
	                                out_form = out_form)
	return(marg_cdf_boot_est)
}

#' rename to reflect that its used for both CDF and PMF now
get_one_marg_cdf <- function(treat, covar, treat_form,
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

  	marg_cdf_est <- marginalize_cdf(cdf_est = cdf_est)
  	marg_pmf_est <- marginalize_pmf(pmf_est = pmf_est)

  	return(list(cdf = marg_cdf_est, pmf = marg_pmf_est))
}

#' this function could be renamed since it's used for CDF and PMF
#' for CDF, set remove_last = FALSE, for PMF set remove_last = TRUE
evaluate_marg_cdf_simul_ci <- function(marg_cdf_est, marg_cdf_eif, alpha, n,
                                       remove_last = FALSE){
	simul_ci <- mapply(pt_est = marg_cdf_est, trt_spec_marg_cdf_eif = marg_cdf_eif, 
	                    FUN = compute_trt_spec_marg_cdf_simul_ci,
	                    SIMPLIFY = FALSE, MoreArgs = list(remove_last = remove_last,
	                                                      alpha = alpha))
    return(simul_ci)
}

#' this function could be renamed since it's used for CDF and PMF
compute_trt_spec_marg_cdf_simul_ci <- function(pt_est, trt_spec_marg_cdf_eif,
                                               remove_last = TRUE, alpha){
	# remove largest value
	if(remove_last){ # for CDF since last pt_est is always 1
		pt_est <- pt_est[-length(pt_est)] 
	}
	K <- length(pt_est)
	gradient <- diag(1 / (pt_est - pt_est^2))
	cor_mat <- cor(trt_spec_marg_cdf_eif %*% gradient)
	# put on logistic scale
	cov_est_logistic <- cov(trt_spec_marg_cdf_eif %*% gradient) / length(trt_spec_marg_cdf_eif[,1])
	# Sigma <- n * cov_est_logistic
	normal_samples <- MASS::mvrnorm(n = 1e5, mu = rep(0, K),
	                                Sigma = cor_mat)
	max_samples <- apply(abs(normal_samples), 1, max)
	q_1alpha <- quantile(max_samples, p = 1 - alpha)
	return(plogis(qlogis(pt_est) + t(c(-q_1alpha, q_1alpha) %o% sqrt(diag(cov_est_logistic)))))
}

evaluate_marg_cdf_ptwise_ci <- function(marg_cdf_est, marg_cdf_eif, alpha){
	marg_cdf_cov <- lapply(marg_cdf_eif, function(x){
		cov(x) / length(x[,1])
	})
	# do on logistic scale
	ptwise_ci <- mapply(pt_est = marg_cdf_est, cov_est = marg_cdf_cov, 
	       				FUN = compute_trt_spec_marg_cdf_ptwise_ci, 
	       				MoreArgs = list(alpha = alpha, cdf = TRUE), SIMPLIFY = FALSE)
    return(ptwise_ci)
}

evaluate_marg_pmf_ptwise_ci <- function(marg_pmf_est, marg_pmf_eif, alpha){
	marg_pmf_cov <- lapply(marg_pmf_eif, function(x){
		cov(x) / length(x[,1])
	})
	# do on logistic scale
	ptwise_ci <- mapply(pt_est = marg_pmf_est, cov_est = marg_pmf_cov, 
	       				FUN = compute_trt_spec_marg_cdf_ptwise_ci, 
	       				MoreArgs = list(alpha = alpha, cdf = FALSE), SIMPLIFY = FALSE)
    return(ptwise_ci)
}

#' this function could be renamed as it's used for both CDF and PMF
compute_trt_spec_marg_cdf_ptwise_ci <- function(pt_est, cov_est, alpha, cdf = TRUE){
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
	 return(plogis(qlogis(pt_est) + t(qnorm(c(alpha/2, 1 - alpha/2)) %o% sqrt(diag(cov_est_logistic)))))
}

evaluate_marg_pmf_eif <- function(pmf_est, treat_prob_est, treat, out, out_levels){
	eif_matrix_list <- mapply(trt_spec_pmf_est = pmf_est, 
	       trt_spec_prob_est = treat_prob_est, trt_level = list(1,0), 
	       FUN = evaluate_trt_spec_pmf_eif,
	       MoreArgs = list(treat = treat, out = out, out_levels = out_levels),
	       SIMPLIFY = FALSE)
	# cov_matrix <- lapply(eif_matrix_list, function(x){
	# 	cov(x) / length(out)
	# })
	return(eif_matrix_list)
}


evaluate_marg_cdf_eif <- function(cdf_est, treat_prob_est, treat, out, out_levels){
	eif_matrix_list <- mapply(trt_spec_cdf_est = cdf_est, 
	       trt_spec_prob_est = treat_prob_est, trt_level = list(1,0), 
	       FUN = evaluate_trt_spec_theta_eif,
	       MoreArgs = list(treat = treat, out = out, out_levels = out_levels),
	       SIMPLIFY = FALSE)
	# cov_matrix <- lapply(eif_matrix_list, function(x){
	# 	cov(x) / length(out)
	# })
	return(eif_matrix_list)
}

marginalize_pmf <- function(pmf_est){
	lapply(pmf_est, colMeans)
}
