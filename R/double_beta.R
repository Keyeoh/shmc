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
#' @return Don't know yet.
#'
#' @importFrom rstan sampling
#' @export
#'
double_beta = function(bs, ox, groups, nu = 1000, alpha = c(1, 1, 1)) {
  dat = list(
    n_probes = nrow(bs),
    n_samples = ncol(bs),
    n_groups = max(groups),
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

  return(fit)
}
