#' @importFrom rlang .data

################################################################################
##### function that returns a parameter based on the log MH ratio supplied #####
################################################################################

acc.ratio <- function(curr_par, prop_par, lr){
  test <- log(stats::runif(1)) < lr
  par_to_return <- curr_par * (1 - test) + prop_par * test
  return(list(par_to_return, as.numeric(test)))
}

################################################################################
## function that checks the acceptance rate and modifies the proposal variance #
################################################################################

change_prop_var <- function(x, iteration){
  acc_rate <- x[1] / iteration
  new_sd <- x[2]
  if(acc_rate < 0.15) new_sd <- x[2] * 0.8
  if(acc_rate >= 0.35 & acc_rate <= 0.70) new_sd <- x[2] * 1.2
  if(acc_rate > 0.70) new_sd <- x[2] * 1.5
  return(new_sd)
}

################################################################################
################# Estimate the Beta Distribution Parameters ####################
################################################################################

init_beta_pars <- function(gene_level_data, binary_ess_ind) {
  mu_var <- plyr::ddply(gene_level_data, ~z_orig, summarise,
                        id_mu = mean(.data$insdens),
                        id_var = stats::var(.data$insdens))
  mu <- mu_var[mu_var$z_orig == binary_ess_ind, 2]
  v <-  mu_var[mu_var$z_orig == binary_ess_ind, 3]
  alpha <- ((1 - mu) / v - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(c(alpha, beta))
}

################################################################################
######################### log likelihood for ESS genes #########################
################################################################################

e_ll <- function(id_e, ae, be) log(stats::dbeta(id_e, ae, be))

################################################################################
######################## log likelihood for NESS genes #########################
################################################################################

n_ll <- function(id_n, an, bn) log(stats::dbeta(id_n, an, bn))

################################################################################
########## calculate p_1 needed to update gene essentiality status #############
################################################################################

calc_logp1 <- function(theta, id, ae, be) log(theta) + e_ll(id, ae, be)

################################################################################
########### calculate p_0 needed to update gene essentiality status ############
################################################################################

calc_logp0 <- function(theta, id, an, bn) log(1 - theta) + n_ll(id, an, bn)
