// =============================================================================
// model_spatial.stan  —  Model 1
// Truncated power-law mixture with spatial weights integrated in the likelihood.
// Spatial weight w_norm encodes proximity to known point sources;
// it is applied only to the astrophysical component.
//
// Parameters:
//   alpha_astro ~ Normal(2.2, 0.5)  — astrophysical spectral index
//   alpha_atmo  ~ Normal(3.7, 0.3)  — atmospheric spectral index
//   pi_astro    ~ Beta(1, 99)       — astrophysical fraction (strong prior toward atm.)
// =============================================================================

data {
  int<lower=1>          N;
  vector[N]             logE;
  vector<lower=0>[N]    w_norm;     // normalised spatial weight per event
  real<lower=0>         E_min;
  real<lower=0>         E_max;
}
transformed data {
  vector[N] E     = pow(10.0, logE);
  vector[N] log_E = log(E);
}
parameters {
  real<lower=1.0>           alpha_astro;
  real<lower=1.0>           alpha_atmo;
  real<lower=0, upper=1>    pi_astro;
}
model {
  // Priors
  alpha_astro ~ normal(2.2, 0.5);
  alpha_atmo  ~ normal(3.7, 0.3);
  pi_astro    ~ beta(1, 99);

  // Likelihood (computed inline; no generated quantities block needed)
  {
    real norm_astro = (alpha_astro - 1.0) /
                      (pow(E_min, 1.0 - alpha_astro) - pow(E_max, 1.0 - alpha_astro));
    real norm_atmo  = (alpha_atmo  - 1.0) /
                      (pow(E_min, 1.0 - alpha_atmo)  - pow(E_max, 1.0 - alpha_atmo));

    vector[N] log_lik_astro = log(w_norm) + log(norm_astro) - alpha_astro * log_E;
    vector[N] log_lik_atmo  =               log(norm_atmo)  - alpha_atmo  * log_E;

    target += sum(log(pi_astro * exp(log_lik_astro) + (1.0 - pi_astro) * exp(log_lik_atmo)));
  }
}
