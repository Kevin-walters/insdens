#' Calculate the posterior probability that a gene is essential using MCMC.
#'
#' @param file_path string containing the path to the source csv file. The
#'  source file must be a csv file with two columns. The first column contains
#'  the gene name and must be called \code{gene}. The second column must contain the
#'  insertion density and must be called \code{insdens}.
#' @param insdens_thresh insertion density threshold used to initialise the
#'  essentiality status of each gene. The default value is 0.025. Use the value
#'  that you think would a priori be your choice of insertion density threshold
#'  to classify genes into essential and non-essential. It must be between 0 and
#'  1.
#' @param init_theta initial value of the theta parameter, the prior probability that a
#'  gene is essential. The default value is 0.05. Use your a priori estimate of the
#'  proportion of essential genes. It must be between 0 and 1.
#' @param niter A positive integer greater than or equal to 2 representing the
#'  total number of MCMC iterations including the \code{burn-in}.
#'  The default value is 21,000.
#' @param burn_in A non-negative integer representing the MCMC \code{burn-in}.The default
#' is 1000. The \code{burn-in} samples are not included in the posterior estimates.
#' It must be less than niter.
#' @param acc_window A positive integer representing the number of steps between
#'  updates to proposal standard deviations. The aim is to achieve an acceptance
#'  rate of between 15\% and 35\%. The default value is 500.
#' @param ae_sd The proposal standard devation for the \eqn{\alpha_E} parameter.
#' The default value is 0.25.
#' @param be_sd The proposal standard devation for the \eqn{\beta_E} parameter.
#' The default value is 0.40.
#' @param an_sd The proposal standard devation for the \eqn{\alpha_N} parameter.
#' The default value is 0.05.
#' @param bn_sd The proposal standard devation for the \eqn{\beta_N} parameter.
#' The default value is 0.05.
#' @param print_prop_sd A logical value indicating whether to print the updated proposal
#'  standard deviations and acceptance rates of each parameter. If TRUE they
#'  will be printed  every \code{acc_window} iterations.
#' @param print_it A logical value indicating whether to print the iteration.
#' number every \code{acc_window} iterations. This is helpful for seeing which
#' iteration the chain is on. For example, if \code{acc_window = 100} and
#' \code{print_it = T} then the iteration number will be printed every 100
#' iterations.
#' @param thin A positive integer that represents the number of samples between
#'  the retained samples. For example \code{thin=10} means retain every tenth
#'  sample. The recommended and default value is 1.
#' @return A list of two elements. The first list element is a data frame in
#'  which the first column gives the gene names and the the second column gives
#'  the posterior probability that the gene is essential. The second list
#'  element is a data frame containing the posterior samples of the model
#'  parameters (including the burn-in). This allows the user to check
#'  convergence of the MCMC chain.
#' @importFrom ttutils isInteger
#' @export

post_prob <- function(file_path = "insertion_density_csv_file",
                      insdens_thresh = 0.025,
                      init_theta = 0.05,
                      niter = 21000,
                      burn_in = 1000,
                      acc_window = 500,
                      ae_sd = 0.25,
                      be_sd = 0.30,
                      an_sd = 0.04,
                      bn_sd = 0.05,
                      print_prop_sd = T,
                      print_it = T,
                      thin = 1){

  ##############################################################################
  ######################### do input parameter checks ##########################
  ##############################################################################

  if(init_theta < 0 | init_theta >  1)
  stop("Your initial estimate of the proportion of essential genes (init_theta)
       must be between 0 and 1\n")

  if(niter <= burn_in)
  stop("The burn in cannot be greater than the number of iterations")

  if(niter < 2 | !ttutils::isInteger(niter))
    stop("niter must be an integer at least 2")
  if(burn_in < 0 | !ttutils::isInteger(burn_in))
    stop("burn-in must be an a non-negative integer")
  if(acc_window < 1 | !ttutils::isInteger(acc_window))
    stop("acc_window must be an integer at at least 1 ")
  if(thin < 1 | !ttutils::isInteger(thin))
    stop("thin must be an integer at least 1")
  if(!is.logical(print_prop_sd)) stop("print_prop_sd can only be TRUE of FALSE")

  ##############################################################################
  ############ read in data and initialise essential status vector #############
  ##############################################################################

  gene_dat <- utils::read.csv(file_path, header = T)
  # check column names
  if(dim(gene_dat)[2] != 2) stop("csv file must have exactly 2 columns")
  if(names(gene_dat)[1] != "gene" | names(gene_dat)[2] != "insdens")
  stop("columns 1 and 2 of the input csv file must be called gene & insdens
       respectively (case sensitive)")
  num_genes <- dim(gene_dat)[1]
  z_orig <- ifelse(gene_dat$insdens < insdens_thresh, 1, 0)
  if(sum(z_orig) == 0 | sum(z_orig) == num_genes)
  stop("Your insertion density threshold (insdens_thresh) either initialises all
  genes as essential or initialises all genes as non-essential.
  Try a different value. \n")
  gene_dat <- cbind(gene_dat, z_orig)

  ##############################################################################
  #### initialize empty paramater vectors and matrices, and fixed parameters ###
  ##############################################################################

  ae_vec <- be_vec <- an_vec <- bn_vec <- theta_vec <- rep(NA, niter)
  z_mat <- matrix(NA, nrow = niter, ncol = num_genes)

  ##############################################################################
  ########## fill in first element/row of parameter vectors/matrices ###########
  ##############################################################################

  e_beta_pars <- init_beta_pars(gene_dat, 1)
  ae <- ae_vec[1] <- ae_shape <- e_beta_pars[1]
  be <- be_vec[1] <- be_shape <- e_beta_pars[2]

  n_beta_pars <- init_beta_pars(gene_dat, 0)
  an <- an_vec[1] <- an_shape <- n_beta_pars[1]
  bn <- bn_vec[1] <- bn_shape <- n_beta_pars[2]
  ae_rate <- be_rate <- an_rate <- bn_rate <- 1
  theta <- theta_vec[1] <- init_theta
  z_vec <- z_mat[1, ] <- z_orig

  accept <- cbind("acc"= rep(0, 4), "sd " = c(ae_sd, be_sd, an_sd, bn_sd))
  #accept <- cbind("acc"= rep(0, 4), "sd " = seq(1,4))
  row.names(accept) <- c("a_e", "b_e", "a_n", "b_n")

  id_e <- gene_dat$insdens[z_orig == 1]
  id_n <- gene_dat$insdens[z_orig == 0]

  ##############################################################################
  ################################# run the mcmc ###############################
  ##############################################################################

  for (i in 2: niter) {
    res1 <- ae_update_MH(ae, be, id_e, ae_shape, ae_rate, prop_sd = accept[1, 2])
    ae_vec[i] <- ae <- res1[[1]]

    res2 <- be_update_MH(ae, be, id_e, be_shape, be_rate, prop_sd = accept[2, 2])
    be_vec[i] <- be <-  res2[[1]]

    res3 <- an_update_MH(an, bn, id_n, an_shape, an_rate, prop_sd = accept[3, 2])
    an_vec[i] <- an <- res3[[1]]

    res4 <- bn_update_MH(an, bn, id_n, bn_shape, bn_rate, prop_sd = accept[4, 2])
    bn_vec[i] <- bn <-  res4[[1]]

    accept[, 1] <- accept[, 1] +  c(res1[[2]], res2[[2]], res3[[2]], res4[[2]])
    theta_vec[i] <- theta <- stats::rbeta(1, sum(z_vec) + 0.1,
                                          sum(1 - z_vec) + 0.9)
    prob_vec <- z_update(theta, gene_dat$insdens, ae, be, an, bn)
    z_mat[i,] <- z_vec <- stats::rbinom (num_genes, size = 1, prob = prob_vec)

    id_e <- gene_dat$insdens[z_vec == 1]
    id_n <- gene_dat$insdens[z_vec == 0]

    if(i%%acc_window == 0){
      if(print_it == T) cat("iteration", i, "\n")
      accept[, 2] <- mapply(FUN = change_prop_var, split(accept, 1:4),
                            iteration = acc_window)
      #cat("proposal sd =", round(accept[, 2], digits = 2),
      #    "\n accept rate=", round(accept[, 1] / i, digit = 2), "\n")
      mat_print <- accept
      mat_print[, 1] <- accept[, 1] / acc_window
      if(print_prop_sd == T) print(mat_print, digits = 2)
      accept[, 1] <- rep(0, 4)
    }
  }
  if(burn_in > 0){
    ae_kept <- ae_vec[-c(1:burn_in)]
    an_kept <- an_vec[-c(1:burn_in)]
    be_kept <- be_vec[-c(1:burn_in)]
    bn_kept <- bn_vec[-c(1:burn_in)]
    theta_kept <- theta_vec[-c(1:burn_in)]
    z_mat_kept <- z_mat[-c(1:burn_in), ]
  } else{
    ae_kept <- ae_vec
    an_kept <- an_vec
    be_kept <- be_vec
    bn_kept <- bn_vec
    theta_kept <- theta_vec
    z_mat_kept <- z_mat
  }

  if(thin > 1){
    its_kept <- niter - burn_in
    ae_kept <- ae_kept[seq(1, its_kept) %% thin == 0]
    an_kept <- an_kept[seq(1, its_kept) %% thin == 0]
    be_kept <- be_kept[seq(1, its_kept) %% thin == 0]
    bn_kept <- bn_kept[seq(1, its_kept) %% thin == 0]
    theta_kept <- theta_kept[seq(1, its_kept) %% thin == 0]
    z_mat_kept <- z_mat_kept[seq(1, its_kept) %% thin == 0, ]
  }

  ppe <- apply(z_mat_kept, 2, mean)
  return (list(data.frame("gene" = gene_dat$gene, "post_prob" = ppe),
               data.frame("alpha_e"= ae_kept, "alpha_n" = an_kept,
                          "beta_e" = be_kept, "beta_n" = bn_kept,
                          "theta" = theta_kept)))
}
