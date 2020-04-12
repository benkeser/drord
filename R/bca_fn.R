bca_interval <- function(pt_est, boot_samples, 
                         jack_samples, alpha = 0.05){
	# the usual normal quantiles
	z_alpha2 <- qnorm(alpha/2)
	z_1alpha2 <- qnorm(1 - alpha/2)
	mean_jack <- mean(jack_samples, na.rm = TRUE)
	z0 <- qnorm(mean(boot_samples < pt_est, na.rm = TRUE))
	a <- sum((jack_samples - mean_jack)^3) / (6 * sum((jack_samples - mean_jack)^2)^(3/2))
	alpha1 <- pnorm(z0 + (z0 + z_alpha2) / (1 - a*(z0 + z_alpha2)))
	alpha2 <- pnorm(z0 + (z0 + z_1alpha2) / (1 - a*(z0 + z_1alpha2)))
	bca_ci <- quantile(boot_samples, p = c(alpha1, alpha2), na.rm = TRUE)
	return(as.numeric(bca_ci))
}

compute_trt_spec_bca_intervals <- function(dist = c("cdf","pmf"),
                                           trt = c(1,0),
                                           marg_est,
                                           boot_samples, 
                                           jack_samples, 
                                           alpha){
	idx <- ifelse(trt == 1, 1, 2)
	boot_mat <- t(apply(boot_samples, 2, function(x){
		unlist(x[[dist]][[idx]], use.names = FALSE)
	}))
	jack_mat <- t(apply(jack_samples, 2, function(x){
		unlist(x[[dist]][[idx]], use.names = FALSE)
	}))
	K <- ncol(boot_mat)
	n_est <- ifelse(dist == "cdf", K-1, K)
	bca_ptwise_ci <- matrix(NA, ncol = 2, nrow = n_est)
	for(k in 1:n_est){
		bca_ptwise_ci[k,] <- plogis(bca_interval(pt_est = qlogis(marg_est[[idx]][k]),
	                            boot_samples = qlogis(boot_mat[,k]),
	                            jack_samples = qlogis(jack_mat[,k]),
	                            alpha = alpha))
	}
	# for simultaneous
	cor_mat <- cor(scale(trimmed_logit(boot_mat[1:n_est, 1:n_est])))
	normal_samples <- MASS::mvrnorm(n = 1e5, mu = rep(0, n_est),
	                                Sigma = cor_mat)
	max_samples <- apply(abs(normal_samples), 1, max)
	q_1alpha <- quantile(max_samples, p = 1 - alpha)
	scale_factor <- q_1alpha / qnorm(1 - alpha/2)
	bca_simul_ci <- cbind(
        plogis(trimmed_logit(marg_est[[idx]][1:n_est]) - (trimmed_logit(marg_est[[idx]][1:n_est]) - trimmed_logit(bca_ptwise_ci[,1]))*scale_factor),
        plogis(trimmed_logit(marg_est[[idx]][1:n_est]) + (trimmed_logit(bca_ptwise_ci[,2]) - trimmed_logit(marg_est[[idx]][1:n_est]))*scale_factor)
    )
    return(list(ptwise = bca_ptwise_ci,
                simul = bca_simul_ci))
		
}
      
trimmed_logit <- function(x){
	if(!is.matrix(x)){
		rslt <- x
		rslt[rslt == 0] <- .Machine$double.neg.eps
		rslt[rslt == 1] <- 1 - .Machine$double.neg.eps
	}else{
		rslt <- apply(x, 2, function(y){
			y[y == 0] <- .Machine$double.neg.eps
			y[y == 1] <- 1 - .Machine$double.neg.eps
			y
		})
	}
	return(qlogis(rslt))
}