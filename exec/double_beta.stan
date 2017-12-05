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

transformed parameters {
  real <lower = 0> BS_alpha[n_groups];
  real <lower = 0> BS_beta[n_groups];
  real <lower = 0> OX_alpha[n_groups];
  real <lower = 0> OX_beta[n_groups];

  for (j in 1:n_groups) {
    OX_alpha[j] = meth_state[j][2] * (nu_minus_one + 1);
    OX_beta[j] = (1 - meth_state[j][2]) * (nu_minus_one + 1);
    BS_alpha[j] = (meth_state[j][2] + meth_state[j][3]) * (nu_minus_one + 1);
    BS_beta[j] = (1 - meth_state[j][2] - meth_state[j][3]) * (nu_minus_one + 1);
  }
}

model {
  nu_minus_one ~ gamma(nu_shape, nu_shape / nu_mean);

  for (j in 1:n_groups) {
    meth_state[j] ~ dirichlet(alpha);
  }

  for (j in 1:n_samples) {
    BS[, j] ~ beta(BS_alpha[group_id[j]], BS_beta[group_id[j]]);
    OX[, j] ~ beta(OX_alpha[group_id[j]], OX_beta[group_id[j]]);
  }
}
