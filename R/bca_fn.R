#' Compute a BCa confidence interval 
#' 
#' @param pt_est The point estimate of the parameter of interest
#' @param boot_samples A collection of bootstrap realizations of the estimator 
#' of the parameter of interest
#' @param jack_samples A vector of jackknife estimates of the 
#' parameter of interest. 
#' @param alpha Confidence intervals have nominal level 1-\code{alpha}. 
#' @return 2-length vector containing BCa confidence interval limits. 
#' @importFrom stats quantile pnorm qnorm
#' 
bca_interval <- function(pt_est, boot_samples, 
                         jack_samples, alpha = 0.05){
	# the usual normal quantiles
	z_alpha2 <- stats::qnorm(alpha/2)
	z_1alpha2 <- stats::qnorm(1 - alpha/2)
	mean_jack <- mean(jack_samples, na.rm = TRUE)
	z0 <- stats::qnorm(mean(boot_samples < pt_est, na.rm = TRUE))
	a <- sum((mean_jack - jack_samples)^3) / (6 * sum((mean_jack - jack_samples)^2)^(3/2))
	alpha1 <- stats::pnorm(z0 + (z0 + z_alpha2) / (1 - a*(z0 + z_alpha2)))
	alpha2 <- stats::pnorm(z0 + (z0 + z_1alpha2) / (1 - a*(z0 + z_1alpha2)))
	bca_ci <- stats::quantile(boot_samples, p = c(alpha1, alpha2), na.rm = TRUE)
	return(as.numeric(bca_ci))
}

#' Used to compute treatment-specific BCa intervals for 
#' the CDF and PMF
#' 
#' @param dist Which one? CDF or PMF?
#' @param trt Which treatment?
#' @param marg_est The point estimate 
#' @param boot_samples A collection of bootstrap realizations of the estimator 
#' of the parameter of interest
#' @param jack_samples A vector of jackknife estimates of the 
#' parameter of interest. 
#' @param alpha Confidence intervals have nominal level 1-\code{alpha}. 
#' @importFrom stats cor
#' @return List of pointwise and simultaneous confidence intervals for \code{dist}.
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
		bca_ptwise_ci[k,] <- stats::plogis(bca_interval(pt_est = stats::qlogis(marg_est[[idx]][k]),
	                            boot_samples = stats::qlogis(boot_mat[,k]),
	                            jack_samples = stats::qlogis(jack_mat[,k]),
	                            alpha = alpha))
	}
	# for simultaneous
	cor_mat <- stats::cor(scale(trimmed_logit(boot_mat[1:n_est, 1:n_est])))
	normal_samples <- MASS::mvrnorm(n = 1e5, mu = rep(0, n_est),
	                                Sigma = cor_mat)
	max_samples <- apply(abs(normal_samples), 1, max)
	q_1alpha <- stats::quantile(max_samples, p = 1 - alpha)
	scale_factor <- q_1alpha / qnorm(1 - alpha/2)
	bca_simul_ci <- cbind(
        stats::plogis(trimmed_logit(marg_est[[idx]][1:n_est]) - (trimmed_logit(marg_est[[idx]][1:n_est]) - trimmed_logit(bca_ptwise_ci[,1]))*scale_factor),
        stats::plogis(trimmed_logit(marg_est[[idx]][1:n_est]) + (trimmed_logit(bca_ptwise_ci[,2]) - trimmed_logit(marg_est[[idx]][1:n_est]))*scale_factor)
    )
    return(list(ptwise = bca_ptwise_ci,
                simul = bca_simul_ci))
		
}
 
#' Trimmed logistic function
#' @param x A numeric between 0 and 1
trimmed_logit <- function(x){
	if(!is.matrix(x)){
		rslt <- x
		rslt[rslt <= 0] <- .Machine$double.neg.eps
		rslt[rslt >= 1] <- 1 - .Machine$double.neg.eps
	}else{
		rslt <- apply(x, 2, function(y){
			y[y <= 0] <- .Machine$double.neg.eps
			y[y >= 1] <- 1 - .Machine$double.neg.eps
			y
		})
	}
	return(stats::qlogis(rslt))
}