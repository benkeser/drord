
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

	rslt <- list(est = wmean_est, cov = NULL)

	if(return_cov){
		eif_mat <- Reduce(cbind, wmean_eif_est)
		colnames(eif_mat) <- NULL
		rslt$cov <- cov(eif_mat) / length(out)
	}
	return(rslt)
}

estimate_cond_mean <- function(trt_spec_pmf_est, ordered_out_levels, ordered_out_weights){
	rslt <- apply(trt_spec_pmf_est, 1, function(pmf){
		sum(ordered_out_levels * ordered_out_weights * pmf) / sum(ordered_out_weights * pmf)
	})
	return(rslt)
}

estimate_eif_wmean <- function(trt_spec_cond_mean_est,
                              trt_spec_prob_est, 
                              trt_level,
                              out, treat){
	eif <- as.numeric(treat == trt_level) / trt_spec_prob_est * (out - trt_spec_cond_mean_est) + 
			trt_spec_cond_mean_est - mean(trt_spec_cond_mean_est)
	return(eif)
}

estimate_ci_wmean <- function(
  out,
  treat,
  covar,
  wmean_est,
  alpha = 0.05, 
  out_levels = order(unique(out)),
  out_form = NULL,
  out_weights = rep(1, length(out_levels)),
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
		                    alpha = alpha, out_weights = out_weights)
	}else{
		bca_ci <- NULL
	}
	return(list(wald = wald_ci, bca = bca_ci))
}

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

#' following slides here: http://users.stat.umn.edu/~helwig/notes/bootci-Notes.pdf
bca_wmean <- function(treat, covar, out, nboot, 
                      treat_form, out_levels, out_form, out_weights,
                      wmean_est, alpha = 0.05){
	boot_samples <- replicate(nboot, 
	                          one_boot_wmean(treat = treat, 
                                             covar = covar, 
	                                         out = out, 
	                                         treat_form = treat_form, 
	                                         out_levels = out_levels, 
	                                         out_form = out_form,
	                                         out_weights = out_weights))
	boot_trt1 <- boot_samples[1,]
	boot_trt0 <- boot_samples[2,]
	boot_diff <- boot_samples[1,] - boot_samples[2,]

	jack_samples <- jack_wmean(treat = treat,
	                           covar = covar,
	                           out = out, 
	                           treat_form = treat_form, 
	                           out_levels = out_levels, 
	                           out_form = out_form,
	                           out_weights = out_weights)

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

jack_wmean <- function(treat, covar, out, treat_form, out_levels, out_form, out_weights){
  	wmean_jack_est <- sapply(seq_along(out), function(i){
		wmean_minusi <- get_one_wmean(treat = treat[-i],
		                                covar = covar[-i, , drop = FALSE],
		                                out = out[-i],
		                                treat_form = treat_form,
		                                out_levels = out_levels,
		                                out_form = out_form,
		                                out_weights = out_weights)  		
		return(wmean_minusi)
  	})
  	return(wmean_jack_est)
}

one_boot_wmean <- function(treat, covar, out, treat_form, out_levels, out_form, out_weights){
	boot_idx <- sample(seq_along(out), replace = TRUE)
	wmean_boot_est <- tryCatch({get_one_wmean(treat = treat[boot_idx],
	                                covar = covar[boot_idx, , drop = FALSE],
	                                out = out[boot_idx],
	                                treat_form = treat_form,
	                                out_levels = out_levels,
	                                out_form = out_form,
	                                out_weights = out_weights)}, error = function(e){
		list(rep(NA, 3))
	})
	return(wmean_boot_est)
}
get_one_wmean <- function(treat, covar, treat_form,
                          out, out_levels, out_form,
                          out_weights){
	# obtain estimate of treatment probabilities
	treat_prob_est <- estimate_treat_prob(treat = treat,
	                                      covar = covar,
	                                      treat_form = treat_form)

	# obtain estimate of conditional PMF for each treatment level
	pmf_est <- estimate_pmf(out = out, treat = treat, 
	                        covar = covar, out_levels = out_levels,
	                        out_form = out_form, treat_prob_est = treat_prob_est)
	
  	wmean_est <- estimate_wmean(
          pmf_est = pmf_est, treat = treat, out = out, out_levels = out_levels, 
          out_weights = out_weights, treat_prob_est = treat_prob_est, 
          return_cov = FALSE
        )
  	return(wmean_est$est)
}