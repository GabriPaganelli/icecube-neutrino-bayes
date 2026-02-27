// =============================================================================
// model_uniform.stan  —  Model 3
// Truncated power-law mixture, energy-only (no spatial weights).
// Uses uniform priors via hard parameter bounds; same log_diff_exp
// normalisation as model_informed.stan.
// Spatial information is combined in a second stage (R/05_posteriors.R).
//
// Parameters:
//   alpha_astro ~ Uniform(1, 4)  — astrophysical spectral index (diffuse prior)
//   alpha_atmo  ~ Uniform(1, 5)  — atmospheric spectral index  (diffuse prior)
//   pi_astro    ~ Beta(2, 3)     — astrophysical fraction (weakly informative)
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
  real<lower=1.0, upper=4.0>    alpha_astro;
  real<lower=1.0, upper=5.0>    alpha_atmo;
  real<lower=0, upper=1>        pi_astro;
}
model {
  // Uniform priors (enforced by parameter bounds above; stated explicitly)
  alpha_astro ~ uniform(1.0, 4.0);
  alpha_atmo  ~ uniform(1.0, 5.0);
  pi_astro    ~ beta(2, 3);

  {
    real log_norm_astro;
    real log_norm_atmo;

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
