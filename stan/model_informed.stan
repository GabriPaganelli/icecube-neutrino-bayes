// =============================================================================
// model_informed.stan  —  Model 2
// Truncated power-law mixture, energy-only (no spatial weights).
// Normalization computed in log space using log_diff_exp for numerical stability.
// Spatial information is combined in a second stage (R/05_posteriors.R).
//
// Parameters:
//   alpha_astro ~ Normal(2.2, 0.5)  — astrophysical spectral index (informed)
//   alpha_atmo  ~ Normal(3.7, 0.3)  — atmospheric spectral index (informed)
//   pi_astro    ~ Beta(2, 35)       — astrophysical fraction (informed, less extreme)
// =============================================================================

data {
  int<lower=1>    N;
  vector[N]       logE;
  real<lower=0>   E_min;
  real<lower=0>   E_max;
}
transformed data {
  vector[N] E       = pow(10.0, logE);
  vector[N] log_E   = log(E);
  real log_E_min    = log(E_min);
  real log_E_max    = log(E_max);
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
  pi_astro    ~ beta(2, 35);

  {
    real log_norm_astro;
    real log_norm_atmo;

    // log-normalisation: log[ (alpha-1) / (E_min^(1-alpha) - E_max^(1-alpha)) ]
    // For alpha > 1: E_min^(1-alpha) > E_max^(1-alpha), so log_diff_exp is valid.
    {
      real exp_term = 1.0 - alpha_astro;
      log_norm_astro = log(alpha_astro - 1.0) -
                       log_diff_exp(exp_term * log_E_min, exp_term * log_E_max);
    }
    {
      real exp_term = 1.0 - alpha_atmo;
      log_norm_atmo  = log(alpha_atmo - 1.0) -
                       log_diff_exp(exp_term * log_E_min, exp_term * log_E_max);
    }

    vector[N] log_lik_astro = log_norm_astro - alpha_astro * log_E;
    vector[N] log_lik_atmo  = log_norm_atmo  - alpha_atmo  * log_E;

    for (n in 1:N) {
      target += log_sum_exp(
        log(pi_astro)       + log_lik_astro[n],
        log(1.0 - pi_astro) + log_lik_atmo[n]
      );
    }
  }
}
