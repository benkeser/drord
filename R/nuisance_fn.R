
#' Estimate probability of receiving each level of treatment
#' Assumes that \code{treat} comes in as a vector of 0's and 1's
#' @return A list where the first element is estimate of Pr(\code{treat} = 1 | \code{covar})
#' for \code{covar} equal to inputted values of \code{covar}
#' and second element is estimate of Pr(\code{treat} = 0 | \code{covar})
#' for \code{covar} equal to inputted values of \code{covar}
#' @importFrom stats as.formula glm binomial
#' 
estimate_treat_prob <- function(treat, covar, treat_form){
	if (length(unique(treat) == 2)){
      fm_treat <- list(stats::glm(
        formula = stats::as.formula(paste0("A~", treat_form)),
        data = data.frame(A = treat, covar), 
        family = stats::binomial()
      ))
      gn_A <- vector(mode = "list", length = 2)
      gn_A[[1]] <- stats::fitted(fm_treat[[1]])
      gn_A[[2]] <- 1 - gn_A[[1]]
    } else {
      stop("multi-level treatments not yet supported")
    } # end multi-level treatment if
    return(gn_A)
}

#' Get a treatment-specific estimate of the conditional PMF. 
#' Essentially this is a wrapper function for \code{fit_trt_spec_reg}. 
#' @importFrom MASS mvrnorm
estimate_pmf <- function(
  out,
  treat,
  covar,
  out_levels,
  out_form = NULL,
  treat_prob_est, 
  ...
  ){
	out <- mapply(trt_level = list(1,0),
	              trt_spec_prob_est = treat_prob_est,
	              FUN = fit_trt_spec_reg,
	              MoreArgs = list(out = out, treat = treat,
	                              covar = covar, out_form = out_form,
	                              out_levels = out_levels),
	              SIMPLIFY = FALSE)
	return(out)
}


#' Helper function to fit a treatment specific outcome regression. 
#' If there are more than 2 observed levels of the outcome for the
#' specified treatment arm, then \code{polr} is used from the \code{MASS}
#' package. Otherwise logistic regression is used. In both cases, 
#' inverse probability of treatment weights are included in the regression. 
#' If there are levels of the outcome that are not observed in this treatment 
#' group, then 0's are added in. The function returns a matrix with named columns
#' corresponding to each outcome (ordered numerically). The entries 
#' represent the estimated covariate-conditional treatment-specific PMF.
#' 
#' @param out 
#' @param treat
#' @param covar
#' @param trt_level
#' @param out_levels
#' @param out_form
#' @param trt_spec_prob_est A vector of estimates of Pr(\code{treat} = \code{trt_level} | \code{covar}).
#' @importFrom MASS polr
fit_trt_spec_reg <- function(
  trt_level,
  trt_spec_prob_est,
  out,
  treat,
  covar,
  out_levels,
  out_form = NULL,
  ...){

	trt_spec_uniq_outcomes <- unique(out[treat == trt_level])
	# only use polr if more than 2 outcome levels observed
	# in treatment arm
  	use_polr <- length(trt_spec_uniq_outcomes) > 2
	
	if(use_polr){
		suppressWarnings( # weights throw unneeded warnings 
			fm_trt <- MASS::polr(
		      formula = stats::as.formula(paste0("factor(out) ~", out_form)),
		      data = data.frame(out = out, 
		                        covar, 
		                        wt = 1/trt_spec_prob_est)[treat == trt_level, , drop = FALSE],
		      weights = wt
		    )
	    )
	 	pmf_treat <- predict(fm_trt, 
	                      newdata = data.frame(out = out, covar,
	                                           wt = 1/trt_spec_prob_est),
	                      type = "probs")
	 	# add in columns of 0's for unobserved outcome levels and re-order things accordingly
	 	obs_out_levels <- colnames(pmf_treat)
	 	not_obs_out_levels <- which(!(out_levels %in% obs_out_levels))
	 	pmf_treat_new <- pmf_treat
	 	if(length(not_obs_out_levels) > 0){
	 		zero_addition <- matrix(0, nrow = length(out), ncol = length(not_obs_out_levels))
	 		pmf_treat_new <- cbind(pmf_treat, zero_addition)
	 		colnames(pmf_treat_new) <- c(colnames(pmf_treat), out_levels[not_obs_out_levels])
	 	}
	 	pmf_treat <- pmf_treat_new[ , order(as.numeric(colnames(pmf_treat_new)))]
	 	row.names(pmf_treat) <- NULL
	}else{
		relabeled_outcome <- as.numeric(out == trt_spec_uniq_outcomes[2])

		suppressWarnings( # weights throw unneeded warnings 
			fm_trt <- stats::glm(
		      formula = stats::as.formula(paste0("factor(out) ~", out_form)),
		      data = data.frame(out = relabeled_outcome, 
		                        covar, 
		                        wt = 1/trt_spec_prob_est)[treat == trt_level, , drop = FALSE],
		      weights = wt, family = binomial()
		    )
	    )

		pred_of_trt_spec_uniq_outcomes2 <- as.numeric(predict(
          fm_trt, type = "response", newdata = data.frame(out = relabeled_outcome,
                                                          covar, wt = 1/trt_spec_prob_est)
        ))
		pmf_treat <- cbind(1 - pred_of_trt_spec_uniq_outcomes2, pred_of_trt_spec_uniq_outcomes2)
		colnames(pmf_treat) <- trt_spec_uniq_outcomes

		# add in columns of 0's for unobserved outcome levels and re-order things accordingly
	 	obs_out_levels <- colnames(pmf_treat)
	 	not_obs_out_levels <- which(!(out_levels %in% obs_out_levels))
	 	pmf_treat_new <- pmf_treat
	 	if(length(not_obs_out_levels) > 0){
	 		zero_addition <- matrix(0, nrow = length(out), ncol = length(not_obs_out_levels))
	 		pmf_treat_new <- cbind(pmf_treat, zero_addition)
	 		colnames(pmf_treat_new) <- c(colnames(pmf_treat), out_levels[not_obs_out_levels])
	 	}
	 	pmf_treat <- pmf_treat_new[ , order(as.numeric(colnames(pmf_treat_new)))]
	}
	return(pmf_treat)
}

#' Map an estimate of the conditional PMF into an estimate of the conditional CDF
estimate_cdf <- function(pmf_est){
	out <- lapply(pmf_est, function(x){
		t(apply(x, 1, cumsum))
	})
	return(out)
}

