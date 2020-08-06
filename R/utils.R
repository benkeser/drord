utils::globalVariables(c("Outcome","Proportion","Intervention",
                         "ptwise_cil","ptwise_ciu","simul_cil",
                         "simul_ciu"))

#' Print the output of a \code{"drord"} object.
#'
#' @param x A \code{"drord"} object
#' @param ci Which confidence interval should be printed. Defaults to BCa,
#' but it BCa was not computed in call to \code{drord}, defaults back to Wald.
#' @param ... Other arguments (not used)
#' @export
#' @method print drord
print.drord <- function(x, ci = "bca", ...) {
  
  rout <- list()
  if("mann_whitney" %in% x$param){
    if("bca" %in% x$ci){
      rout$mann_whitney <- c(x$mann_whitney$est, x$mann_whitney$ci$bca)
      names(rout$mann_whitney) <- c("est", "bca_cil", "bca_ciu")
    }else if("wald" %in% x$ci){
      rout$mann_whitney <- c(x$mann_whitney$est, x$mann_whitney$ci$wald) 
      names(rout$mann_whitney) <- c("est", "wald_cil", "wald_ciu")
    }else{
      rout$mann_whitney <- c(x$mann_whitney$est, NA, NA)
      names(rout$mann_whitney) <- c("est", "cil", "ciu")
    }
  }

  if("log_odds" %in% x$param){
    if("bca" %in% x$ci){
      rout$log_odds <- data.frame(est = x$log_odds$est, 
                                  bca_cil = x$log_odds$ci$bca[,1],
                                  bca_ciu = x$log_odds$ci$bca[,2])
    }else if("wald" %in% x$ci){
      rout$log_odds <- data.frame(est = x$log_odds$est, 
                                  wald_cil = x$log_odds$ci$wald[,1],
                                  wald_ciu = x$log_odds$ci$wald[,2])      
    }else{
      rout$log_odds <- data.frame(est = x$log_odds$est, 
                                  cil = rep(NA, 3),
                                  ciu = rep(NA, 3))
    }
    row.names(rout$log_odds) <- c("treat1", "treat0", "diff")
  }

  if("weighted_mean" %in% x$param){
    if("bca" %in% x$ci){
      rout$weighted_mean <- data.frame(est = x$weighted_mean$est$est, 
                                  bca_cil = x$weighted_mean$ci$bca[,1],
                                  bca_ciu = x$weighted_mean$ci$bca[,2])
    }else if("wald" %in% x$ci){
      rout$weighted_mean <- data.frame(est = x$weighted_mean$est$est, 
                                  wald_cil = x$weighted_mean$ci$wald[,1],
                                  wald_ciu = x$weighted_mean$ci$wald[,2])      
    }else{
      rout$weighted_mean <- data.frame(est = x$weighted_mean$est$est, 
                                  cil = rep(NA, 3),
                                  ciu = rep(NA, 3))
    }
    row.names(rout$weighted_mean) <- c("treat1", "treat0", "diff")
  }

  print(rout)
  return(invisible(rout))
}

#' Print the output of a \code{"drord"} object.
#'
#' @param x A \code{"drord"} object.
#' @param treat_labels Labels for the treatment variables (treat = 1 followed by treat = 0).
#' @param out_labels Labels for the ordered outcome levels. If \code{dist = "cdf"}, the
#' highest level of outcome will be dropped. 
#' @param dist Which distribution to plot. Valid options are \code{"cdf"} or \code{"pmf"}.
#' @param ... Other arguments (not used)
#' @export
#' @return A list with named entries \code{plot} (a \code{ggplot2} object) and \code{plot_data},
#' the \code{data.frame} from which the plot is made. The latter is included for additional
#' modifications to the plot that are desired. 
#' @importFrom ggplot2 ggplot aes geom_bar scale_y_continuous geom_errorbar 
#' @importFrom ggplot2 theme_bw ggtitle position_dodge
#' @method plot drord
plot.drord <- function(x, 
                       treat_labels = c(1, 0),                        
                       dist = "pmf",
                       out_labels = if(dist == "pmf"){
                                          x$out_levels 
                                        }else{
                                          x$out_levels[-length(x$out_levels)]
                                        }, 
                       ...){
  stopifnot(x$est_dist)
  if(dist == "pmf"){
    stopifnot(length(out_labels) == length(x$out_levels))
  }else{
    stopifnot(length(out_labels) == length(x$out_levels) - 1)
  }
  stopifnot(length(treat_labels) == 2)

  n_grp <- length(out_labels)
  suppressWarnings(
    plot_data <- data.frame(
      Proportion = c(x[[dist]]$est[[1]][1:n_grp], x[[dist]]$est[[2]][1:n_grp]),
      Intervention = c(rep(treat_labels[1], n_grp), rep(treat_labels[2], n_grp)), 
      Outcome = rep(out_labels, 2),
      ptwise_cil = c(x[[dist]]$ci$wald[[1]]$ptwise[,1], x[[dist]]$ci$wald[[2]]$ptwise[,1]),
      ptwise_ciu = c(x[[dist]]$ci$wald[[1]]$ptwise[,2], x[[dist]]$ci$wald[[2]]$ptwise[,2]),
      simul_cil = c(x[[dist]]$ci$wald[[1]]$simul[,1], x[[dist]]$ci$wald[[2]]$simul[,1]),
      simul_ciu = c(x[[dist]]$ci$wald[[1]]$simul[,2], x[[dist]]$ci$wald[[2]]$simul[,2])
    )
  )

  plot_data$Outcome <- factor(plot_data$Outcome, 
                              levels = out_labels, ordered = TRUE)
  plot_data$Intervention <- factor(plot_data$Intervention)
  bar_plot <- ggplot2::ggplot(plot_data, 
                              ggplot2::aes(x = Outcome, y = Proportion, 
                                           fill = Intervention)) + 
              ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge()) + 
              ggplot2::scale_y_continuous(limits = c(0,1)) + 
              ggplot2::geom_errorbar(position = ggplot2::position_dodge(0.8), 
                                     ggplot2::aes(ymin = ptwise_cil, ymax = ptwise_ciu),
                                     width = 0.125) + 
              ggplot2::geom_errorbar(position = ggplot2::position_dodge(1.1), 
                                     aes(ymin = simul_cil, ymax = simul_ciu),
                                     width = 0.125, colour = "gray50") + 
              ggplot2::theme_bw() 
  return(list(plot = bar_plot, plot_data = plot_data))
}


#' Simulated COVID-19 outcomes for hospitalized patients. 
#'
#' A simulated dataset containing outcomes, (hypothetical) treatment,
#' and age group 
#'
#' @format A data frame with 500 rows and 3 variables:
#' \describe{
#'   \item{out}{study outcome, here 1 represents death, 2 intubation, 3 no adverse outcome}
#'   \item{age_grp}{age category with 1 the youngest and 7 the oldest}
#'   \item{treat}{hypothetical treatment, here 1 represents an (effective) active treatment and 0 a control}
#' }
"covid19"