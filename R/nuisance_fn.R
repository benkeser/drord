#' Estimate probability of receiving each level of treatment
#' @param treat A \code{numeric} vector containing treatment status. Only values of
#' 0 or 1 are treated as actual treatment levels. Any other value is assumed to encode
#' a value for which the outcome is missing. 
#' @param covar A \code{data.frame} containing the covariates to include in the working
#' proportional odds model. 
#' @param treat_form The right-hand side of a regression formula for the working model of
#' treatment probability as a function of covariates
#' @param return_models If \code{TRUE} the fitted working proportional odds models
#' and treatment probability models are returned. 
#' @return A list where the first element is estimate of Pr(\code{treat} = 1 | \code{covar})
#' for \code{covar} equal to inputted values of \code{covar}
#' and second element is estimate of Pr(\code{treat} = 0 | \code{covar})
#' for \code{covar} equal to inputted values of \code{covar}
#' @importFrom stats as.formula glm binomial
#' 
estimate_treat_prob <- function(treat, covar, treat_form, return_models){
	fm_treat1 <- stats::glm(
	  formula = stats::as.formula(paste0("treat1 ~", treat_form)),
	  data = data.frame(treat1 = as.numeric(treat == 1), covar), 
	  family = stats::binomial()
	)
	gn_A <- vector(mode = "list", length = 2)
	gn_A[[1]] <- predict(fm_treat1, newdata = data.frame(treat = treat, covar),
	                     type = "response")
	if(all(treat %in% c(0,1))){
		gn_A[[2]] <- 1 - gn_A[[1]]
	}else{
		fm_treat0 <- stats::glm(
	  	  formula = stats::as.formula(paste0("treat0 ~", treat_form)),
	  	  data = data.frame(treat0 = as.numeric(treat == 0), covar), 
		  family = stats::binomial()
		)
		gn_A[[2]] <- predict(fm_treat0, newdata = data.frame(treat = treat, covar),
	                     type = "response")
	}
	rout <- list(gn = gn_A,
                 fm = list(treat1 = fm_treat1,
                           treat0 = NULL))
	if(!all(treat %in% c(0,1))){
		rout$fm$treat0 = fm_treat0
	}

    return(rout)
}

#' Get a treatment-specific estimate of the conditional PMF. 
#' Essentially this is a wrapper function for \code{fit_trt_spec_reg}, which
#' fits the proportion odds model in a given treatment arm. 
#' @param out A \code{numeric} vector containing the outcomes. Missing outcomes are 
#' allowed. 
#' @param treat A \code{numeric} vector containing treatment status. Missing
#' values are not allowed unless the corresponding entry in \code{out} is also missing. 
#' Only values of 0 or 1 are treated as actual treatment levels. Any other value is assumed 
#' to encode a value for which the outcome is missing and the corresponding outcome value is 
#' ignored. 
#' @param covar A \code{data.frame} containing the covariates to include in the working
#' proportional odds model. 
#' @param out_levels A \code{numeric} vector containing all ordered levels of the 
#' outcome. 
#' @param out_form The right-hand side of a regression formula for the working proportional 
#' odds model. NOTE: THIS FORMULA MUST NOT SUPPRESS THE INTERCEPT. 
#' @param out_model Which R function should be used to fit the proportional odds 
#' model. Options are \code{"polr"} (from the \code{MASS} package), 
#' "vglm" (from the \code{VGAM} package), or \code{"clm"} (from the \code{ordinal} package).
#' @param treat_prob_est Estimated probability of treatments, output from call
#' to \code{estimate_treat_prob}.
#' @param return_models If \code{TRUE} the fitted working proportional odds models
#' and treatment probability models are returned. 
#' @param stratify Boolean indicating whether to use nonparametric maximum likelihood
#' (i.e., a stratified estimator). If \code{out_form = "1"}, then a covariate-unadjusted
#' estimate is computed. 
#' @param ... Other options (not used). 
#' @return A list with \code{fm} the fitted model for treatment 1 and 0 (or, if 
#' \code{!return_models} then \code{NULL}) and \code{pmf} the estimated PMF 
#' under treatment 1 and 0 evaluated on each observation. 
#' @importFrom MASS mvrnorm
#' @importFrom VGAM propodds
#' @importFrom stats predict qlogis
estimate_pmf <- function(
  out,
  treat,
  covar,
  out_levels,
  out_form = NULL,
  out_model,
  treat_prob_est, 
  stratify = FALSE,
  return_models = TRUE,
  ...
  ){
	out <- mapply(trt_level = list(1,0),
	              trt_spec_prob_est = treat_prob_est,
	              FUN = fit_trt_spec_reg,
	              MoreArgs = list(out = out, treat = treat,
	                              out_model = out_model[1],
	                              covar = covar, out_form = out_form,
	                              out_levels = out_levels,
	                              stratify = stratify),
	              SIMPLIFY = FALSE)
	rout <- list(fm = list(treat1 = out[[1]]$fm_treat, treat0 = out[[2]]$fm_treat),
	             pmf = list(out[[1]]$pmf_treat, out[[2]]$pmf_treat)) 
	return(rout)
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
#' @param out A \code{numeric} vector containing the outcomes. Missing outcomes are 
#' allowed. 
#' @param treat A \code{numeric} vector containing treatment status. Missing
#' values are not allowed unless the corresponding entry in \code{out} is also missing. 
#' Only values of 0 or 1 are treated as actual treatment levels. Any other value is assumed 
#' to encode a value for which the outcome is missing and the corresponding outcome value is 
#' ignored. 
#' @param covar A \code{data.frame} containing the covariates to include in the working
#' proportional odds model. 
#' @param trt_level Which level of treatment to fit the proportional odds model for
#' @param out_levels A \code{numeric} vector containing all ordered levels of the 
#' outcome. 
#' @param out_form The right-hand side of a regression formula for the working proportional 
#' odds model. NOTE: THIS FORMULA MUST NOT SUPPRESS THE INTERCEPT. 
#' @param out_model Which R function should be used to fit the proportional odds 
#' model. Options are \code{"polr"} (from the \code{MASS} package), 
#' "vglm" (from the \code{VGAM} package), or \code{"clm"} (from the \code{ordinal} package).
#' @param trt_spec_prob_est A vector of estimates of Pr(\code{treat} = \code{trt_level} | \code{covar}).
#' @param stratify Boolean indicating whether to use nonparametric maximum likelihood
#' (i.e., a stratified estimator). If \code{out_form = "1"}, then a covariate-unadjusted
#' estimate is computed. 
#' @param ... Other options (not used).
#' @importFrom MASS polr
#' @importFrom stats model.matrix
#' @importFrom VGAM vglm
#' @importFrom ordinal clm


fit_trt_spec_reg <- function(
  trt_level,
  trt_spec_prob_est,
  out,
  treat,
  covar,
  out_levels,
  out_form = NULL,
  out_model,
  stratify, 
  ...){
	if(!stratify){
		trt_spec_uniq_outcomes <- unique(out[treat == trt_level])
		# only use polr if more than 2 outcome levels observed
		# in treatment arm
	  	multi_level <- length(trt_spec_uniq_outcomes) > 2
		
		if(multi_level){
			if(out_model == "polr"){
				suppressWarnings( # weights throw unneeded warnings 
					fm_trt <- tryCatch({MASS::polr(
				      formula = stats::as.formula(paste0("factor(out, ordered = TRUE) ~", out_form)),
				      data = data.frame(out = out, covar),
				      weights = as.numeric(treat == trt_level)/trt_spec_prob_est
				    )}, error = function(e){
				    	mod_mat <- stats::model.matrix(stats::as.formula(paste0("factor(out) ~", out_form)),
				    	                        data = data.frame(out = out, covar))
				    	n_par <- (dim(mod_mat)[2] - 1) + (length(out_levels) - 1)
						MASS::polr(
					      formula = stats::as.formula(paste0("factor(out) ~", out_form)),
					      data = data.frame(out = out, covar),				      	  
					      weights = as.numeric(treat == trt_level)/trt_spec_prob_est,
					      start = c(rep(-1, length(out_levels) - 1), rep(0, dim(mod_mat)[2]-1))
				    	)
				    })
			    )
			 	pmf_treat <- stats::predict(fm_trt, 
			                      newdata = data.frame(out = out, covar,
			                                           wt = 1/trt_spec_prob_est),
			                      type = "probs")
		 	}else if(out_model == "vglm"){
		 		fm_trt <- VGAM::vglm(formula = stats::as.formula(paste0("factor(out, ordered = TRUE) ~", out_form)),
		 		               family = VGAM::propodds,
		 		               data = data.frame(out = out, covar)[treat == trt_level, , drop = FALSE],
				      		   weights = 1/trt_spec_prob_est[treat == trt_level])
		 		# re-order things
			 	tmp <- stats::predict(fm_trt, newdata = covar, type = "response")
				colnames(tmp) <- sort(trt_spec_uniq_outcomes)[-1]
				pmf_treat <- t(apply(tmp, 1, function(x){
					rr <- rev(diff(c(0,rev(stats::plogis(x)))))
					c(1 - sum(rr), rr)
				}))
				colnames(pmf_treat)[1] <- sort(trt_spec_uniq_outcomes)[1]
		 	}else if(out_model == "clm"){
		 		out_f <- factor(out, ordered = TRUE)
		 		fm_trt <- ordinal::clm(formula = stats::as.formula(paste0("out ~", out_form)),		 		               
		 		               data = data.frame(out = out_f, covar)[treat == trt_level, , drop = FALSE],
				      		   weights = 1/trt_spec_prob_est[treat == trt_level])
			 	pmf_treat <- stats::predict(fm_trt, newdata = data.frame(covar),
		 			        		 type = "prob")$fit
		 	}else if(out_model == "pooled-logistic"){
		 		out_f <- factor(out, ordered = TRUE)
		 		fm_trt <- POplugin(form = stats::as.formula(paste0("factor(out, ordered = TRUE) ~", out_form)),
		 		                   data = data.frame(out = out_f, covar)[treat == trt_level, , drop = FALSE], 
		 		                   weights = 1 / trt_spec_prob_est[treat == trt_level])
		 		pmf_treat <- stats::predict(fm_trt, newdata = data.frame(covar))
		 	}
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
			                        covar)[treat == trt_level, , drop = FALSE],
			      weights = 1/trt_spec_prob_est[treat == trt_level], family = binomial()
			    )
		    )

			pred_of_trt_spec_uniq_outcomes2 <- as.numeric(stats::predict(
	          fm_trt, type = "response", newdata = data.frame(out = relabeled_outcome,
	                                                          covar)
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
	}else{
		fm_trt <- NULL
		if(out_form != "1"){
			if(dim(covar)[2] > 1){
				stop("stratified estimator only implemented for single covariate")
			}
			covar_levels <- unique(covar[,1])
			pmf_treat <- matrix(NA, ncol = length(out_levels), nrow = length(out))
			ct <- 0
			for(o_lev in out_levels[-length(out_levels)]){
				ct <- ct + 1
				for(covar_lev in covar_levels){
					# check that some have this level
					obs_level <- sum(covar[,1] == covar_lev)
					if(obs_level > 0){
						pmf_treat[covar[,1] == covar_lev, ct] <- mean(out[treat == trt_level & covar[,1] == covar_lev] == o_lev)
					}else{
						# if no one observed, just replace with global mean
						# would be better probably to put a kernel thing here... for later.
						pmf_treat[covar[,1] == covar_lev, ct] <- mean(out[treat == trt_level] == o_lev)
					}
				}
			}
			pmf_treat[,ncol(pmf_treat)] <- 1 - rowSums(pmf_treat[,1:(ncol(pmf_treat) - 1)])
		}else{
			if(dim(covar)[2] > 1){
				stop("stratified estimator only implemented for single covariate")
			}
			covar_levels <- unique(covar[,1])
			pmf_treat <- matrix(NA, ncol = length(out_levels), nrow = length(out))
			ct <- 0
			for(o_lev in out_levels[-length(out_levels)]){
				ct <- ct + 1
				for(covar_lev in covar_levels){
					# check that some have this level
					# !!! lazy coding
					obs_level <- 0
					if(obs_level > 0){
						pmf_treat[covar[,1] == covar_lev, ct] <- mean(out[treat == trt_level & covar[,1] == covar_lev] == o_lev)
					}else{
						# if no one observed, just replace with global mean
						# would be better probably to put a kernel thing here... for later.
						pmf_treat[covar[,1] == covar_lev, ct] <- mean(out[treat == trt_level] == o_lev)
					}
				}
			}
			pmf_treat[,ncol(pmf_treat)] <- 1 - rowSums(pmf_treat[,1:(ncol(pmf_treat) - 1)])
		}
	}
	return(list(pmf_treat = pmf_treat,
	            fm_treat = fm_trt))
}

#' Map an estimate of the conditional PMF into an estimate of the conditional CDF
#' @param pmf_est A list of the treatment-specific PMF estimates
#' @return A list of treatment-specific CDF estimates
estimate_cdf <- function(pmf_est){
	out <- lapply(pmf_est, function(x){
		t(apply(x, 1, cumsum))
	})
	return(out)
}
