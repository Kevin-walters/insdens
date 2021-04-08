################################################################################
######################### function for updating ae ##########################
################################################################################

ae_update_MH <- function(ae, be, id_e, ae_shape, ae_rate, prop_sd) {
  prop_ae <- exp(stats::rnorm(1, 0, prop_sd)) * ae
  curr_ess_ll <- sum(e_ll(id_e, ae, be))
  prop_ess_ll <- sum(e_ll(id_e, prop_ae, be))
  prop_log_cond <- log(stats::dgamma(prop_ae, ae_shape, ae_rate)) + prop_ess_ll
  curr_log_cond <- log(stats::dgamma(ae, ae_shape, ae_rate)) + curr_ess_ll
  lr <- prop_log_cond - curr_log_cond
  acc.ratio(ae, prop_ae, lr)
}

################################################################################
########################## function for updating be #########################
################################################################################

be_update_MH <- function(ae, be, id_e, be_shape, be_rate, prop_sd) {
  prop_be <- exp(stats::rnorm(1, 0, prop_sd)) * be
  curr_ll <- sum(e_ll(id_e, ae, be))
  prop_ll <- sum(e_ll(id_e, ae, prop_be))
  prop_log_cond <- log(stats::dgamma(prop_be, be_shape, be_rate)) + prop_ll
  curr_log_cond <- log(stats::dgamma(be, be_shape, be_rate)) + curr_ll
  lr <- prop_log_cond - curr_log_cond
  acc.ratio(be, prop_be, lr)
}

################################################################################
######################### function for updating an #########################
################################################################################

an_update_MH <- function(an, bn, id_n, an_shape, an_rate, prop_sd) {
  prop_an <- exp(stats::rnorm(1, 0, prop_sd)) * an
  curr_ll <- sum(n_ll(id_n, an, bn))
  prop_ll <- sum(n_ll(id_n, prop_an, bn))
  prop_log_cond <- log(stats::dgamma(prop_an, an_shape, an_rate)) + prop_ll
  curr_log_cond <- log(stats::dgamma(an, an_shape, an_rate)) + curr_ll
  lr <- prop_log_cond - curr_log_cond
  acc.ratio(an, prop_an, lr)
}

################################################################################
########################## function for updating bn ########################
################################################################################

bn_update_MH <- function(an, bn, id_n, bn_shape, bn_rate, prop_sd) {
  prop_bn <- exp(stats::rnorm(1, 0, prop_sd)) * bn
  curr_ll <- sum(n_ll(id_n, an, bn))
  prop_ll <- sum(n_ll(id_n, an, prop_bn))
  prop_log_cond <- log(stats::dgamma(prop_bn, bn_shape, bn_rate)) + prop_ll
  curr_log_cond <- log(stats::dgamma(bn, bn_shape, bn_rate)) + curr_ll
  lr <- prop_log_cond - curr_log_cond
  acc.ratio(bn, prop_bn, lr)
}

################################################################################
############################ function for updating z ###########################
################################################################################

z_update <- function(theta, id, ae, be, an, bn){
  log_p1 <- calc_logp1(theta, id, ae, be)
  log_p0 <- calc_logp0(theta, id, an, bn)
  1/(1 + exp(log_p0 - log_p1))
}
