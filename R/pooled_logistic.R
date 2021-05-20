globalVariables("wgts")

#' Get a response from model formula
#' 
#' @param formula The model formula
#' @param data The data frame associated with the model
#' @importFrom stats terms as.formula
getResponseFromFormula = function(formula, data) {
    if (attr(stats::terms(stats::as.formula(formula), data = data), 
             which = 'response'))
        return(all.vars(formula)[1])
    else
        return(NULL)
}

#' Fits a proportional odds model via pooled logistic regression.
#' 
#' The outcome in \code{data} (indicated in the \code{form} object) should be an ordered factor.
#' 
#' @importFrom stats model.extract model.frame glm
#' @param form The model formula
#' @param data The data set used to fit the model 
#' @param weights Either equal to 1 (no weights) or a vector of length equal to nrow(data)
#' @return A list with the fitted glm, the original data, levels of the outcome, and the outcome name
#' @importFrom stats update formula
POplugin = function(form, data, weights = 1){
	# name of outcome from the provided formula
	out_name = getResponseFromFormula(form, data) 

	# name of new (binary) outcome column
	new_out_name = paste0(out_name,"Binary") 

	# make sure that this new column name doesn't already appear in the dataframe. If it does, append "1" to it until it doesn't anymore
	while(new_out_name %in% colnames(data)){ 
		new_out_name = paste0(new_out_name,"1")
	}

	outcomes = stats::model.extract(stats::model.frame(form,data),"response") # extract column with outcomes from data

	data_pooled = do.call(rbind,lapply(levels(outcomes),function(lev){ # create pooled dataset to be used when fitting the model
		if(lev >= min(outcomes) & lev < max(outcomes)){
			out = data.frame(
				outcomes<=lev,
				data[,colnames(data)!=out_name,drop=FALSE],
				wgts = weights,
				lev)
			colnames(out)[c(1,ncol(out))] = c(new_out_name,out_name)
			return(out) 
		}
	}))

	glm_fit = suppressWarnings({
		stats::glm(stats::update(form,stats::formula(paste0(new_out_name,"~.-1+",out_name))), data = data_pooled, family=binomial(),weights=wgts)
	})

	out = list(glm_fit=glm_fit,data=data,levs=levels(outcomes),out_name=out_name)
	class(out) <- "POplugin"

	return(out)
}

#' Predict method for a \code{POplugin} object
#' 
#' @param object An object of class \code{POplugin}
#' @param newdata A \code{data.frame} on which to predict
#' @return A data frame with nrow = number of rows in newdata (or the orignal data frame)
#' and with the number of columns equal to the number of levels of the outcome observed in the original data frame
predict.POplugin = function(object,newdata=NULL){
	if(is.null(newdata)){
		newdata = object$data
	}
	levs = object$levs
	out_name = object$out_name
	cumulative = sapply(levs[-length(levs)],function(lev){	
		df = data.frame(
				newdata[,colnames(newdata)!=out_name,drop=FALSE],
				lev)
		colnames(df)[ncol(df)] = out_name
		stats::predict(object$glm_fit,newdata=df,type="response")
		})
	pmf = cbind(cumulative,1) - cbind(0,cumulative)
	colnames(pmf) = levs

	return(pmf)
}

