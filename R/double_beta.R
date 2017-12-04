#' Fit a double-beta model to 5hmC data
#'
#' This function fits a double-beta model to the 5hmC data represented with
#' two matrices containing the Bisulphite and Oxidated data respectively. This
#' model treats nu and alpha as user-fixed parameters.
#'
#' @param bs A p x n numerical matrix containing the Bisulphite data.
#' @param ox A p x n numerical matrix containing the Oxidated data.
#' @param groups A vector of dimension n representing the correspondence between
#' samples (columns) and groups.
#' @param nu A numeric estimate for the total intensity of an average probe.
#' @param alpha A vector of dimension 3 containing the prior estimate for the
#' proportions between no-mc, 5mC and 5hmC.
#'
#' @return A data.frame containing the tidy results from the model fit.
#'
#' @importFrom dplyr %>% select mutate bind_rows everything
#' @importFrom rstan sampling
#' @export
#'
double_beta = function(bs, ox, groups, nu = 1000, alpha = c(1, 1, 1)) {

  n_groups = max(groups)

  dat = list(
    n_probes = nrow(bs),
    n_samples = ncol(bs),
    n_groups = n_groups,
    group_id = groups,
    BS = bs,
    OX = ox,
    nu = nu,
    alpha = alpha
  )

  fit = sampling(
    object = stanmodels$double_beta,
    data = dat
  )

  mc_result = resume_double_beta(fit, n_groups, mod = 'mc')
  hmc_result = resume_double_beta(fit, n_groups, mod = 'hmc')

  result = bind_rows(
    mc_result %>%
      mutate(Mod = '5mC') %>%
      select(Mod, everything()),
    hmc_result %>%
      mutate(Mod = '5hmC') %>%
      select(Mod, everything())
  )

  return(result)
}

#' Resumes a double beta model.
#'
#' This function takes a fitted double beta model as input and generates a tidy
#' data.frame containing all the relevant information from the sampling process.
#'
#' @param fitted_model A stanfit object containing the fitted model.
#' @param n_grp A numerical containing the number of groups in the model.
#' @param mod A character indicating the modification we want to resume.
#'
#' @return A tidy data.frame containing the results of the sampling.
#'
#' @importFrom coda effectiveSize HPDinterval
#' @importFrom dplyr %>% group_by summarize mutate
#' @importFrom tidyr gather
#' @importFrom rstan summary
#' @importFrom stats sd setNames
#' @keywords internal
#'
resume_double_beta = function(fitted_model, n_grp, mod = c('mc', 'hmc')) {
  mod = match.arg(mod)
  mod_id = ifelse(mod == 'mc', '2', '3')
  param_names = paste0('meth_state[', 1:n_grp, ',', mod_id, ']')

  mc_draws = fitted_model %>%
    as.data.frame(pars = param_names) %>%
    setNames(paste0('Group_', 1:n_grp))

  mc_dif_draws = as.matrix(mc_draws) %*% pair_contrasts(n_grp)
  mc_all_draws = cbind(mc_draws, mc_dif_draws)

  mc_rhat = setNames(
    summary(
      fitted_model,
      pars = param_names
    )$summary[, 'Rhat'],
    paste0('Group_', 1:n_grp)
  )

  mc_n_eff = setNames(
    summary(
      fitted_model,
      pars = param_names
    )$summary[, 'n_eff'],
    paste0('Group_', 1:n_grp)
  )

  mc_result = mc_all_draws %>%
    gather(Contrast, Value) %>%
    group_by(Contrast) %>%
    summarize(
      Mean = mean(Value),
      Sd = sd(Value),
      ESS = effectiveSize(Value),
      Lower = HPDinterval(coda::as.mcmc(Value))[, 'lower'],
      Upper = HPDinterval(coda::as.mcmc(Value))[, 'upper']
    ) %>%
    mutate(Std_Err = Sd / sqrt(ESS)) %>%
    mutate(Rhat = mc_rhat[Contrast]) %>%
    mutate(N_Eff = mc_n_eff[Contrast])
}

#' Generates a contrasts matrix
#'
#' This function generates a matrix containing the contrasts corresponding to
#' all the possible pairwise comparisons between the number of groups given as
#' input.
#'
#' @param n A numerical indicating the number of possible groups.
#'
#' @return A matrix containing the pairwise contrasts.
#'
#' @importFrom utils combn
#' @keywords internal
#'
pair_contrasts = function(n) {
  combs = combn(1:n, 2)
  result = matrix(0, nrow = n, ncol = ncol(combs))
  cnames = rep('', ncol = ncol(combs))

  for (i in 1:ncol(combs)) {
    from_idx = combs[1, i]
    result[from_idx, i] = -1

    to_idx = combs[2, i]
    result[to_idx, i] = +1

    cnames[i] = paste0('Group_', to_idx, ' - Group_', from_idx)
  }
  colnames(result) = cnames
  return(result)
}
