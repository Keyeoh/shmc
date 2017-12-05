data {
  int <lower = 1> n_probes;
  int <lower = 1> n_samples;
  int <lower = 1> n_groups;
  real <lower = 0, upper = 1> BS[n_probes, n_samples];
  real <lower = 0, upper = 1> OX[n_probes, n_samples];
  int <lower = 1> group_id[n_samples];
  vector <lower = 1>[3] alpha;
  real <lower = 0> nu_shape;
  real <lower = 0> nu_mean;
}

parameters {
  real <lower = 0> nu_minus_one;
  simplex[3] meth_state[n_groups];
}

model {
  vector[n_groups] BS_alpha;
  vector[n_groups] BS_beta;
  vector[n_groups] OX_alpha;
  vector[n_groups] OX_beta;
  vector[n_groups] mc;
  vector[n_groups] hmc;
  real nu;

  nu_minus_one ~ gamma(nu_shape, nu_shape / nu_mean);
  nu = nu_minus_one + 1;

  for (j in 1:n_groups) {
    meth_state[j] ~ dirichlet(alpha);

    mc[j] = meth_state[j][2];
    hmc[j] = meth_state[j][3];
  }

  OX_alpha = mc * nu;
  OX_beta = (1 - mc) * nu;
  BS_alpha = (mc + hmc) * nu;
  BS_beta = (1 - mc - hmc) * nu;

  for (j in 1:n_samples) {
    BS[, j] ~ beta(BS_alpha[group_id[j]], BS_beta[group_id[j]]);
    OX[, j] ~ beta(OX_alpha[group_id[j]], OX_beta[group_id[j]]);
  }
}
