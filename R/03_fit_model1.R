# =============================================================================
# 03_fit_model1.R
# PURPOSE : Fit Model 1 — truncated power-law energy model with spatial
#           weights integrated directly into the Stan likelihood.
#           Priors: alpha_astro ~ N(2.2, 0.5), alpha_atmo ~ N(3.7, 0.3),
#                   pi_astro ~ Beta(1, 99).
# INPUT   : output/data.rds, stan/model_spatial.stan
# OUTPUT  : output/fit_model1.rds  — Stan fit object
#           output/gamma_model1.rds — posterior classification probabilities
# RUNTIME : ~ 8 hours  (4 chains × 4500 iter, 2 cores)
# WARNING : This script triggers a long MCMC run. Pre-computed results are
#           already in output/. Run only if you want to reproduce the fit.
# =============================================================================

library(rstan)

options(mc.cores = 2)
rstan_options(auto_write = TRUE)

SAVE_BACKUP <- FALSE   # set to TRUE to write a timestamped backup of the fit

# -- Prepare data -------------------------------------------------------------

data_full <- readRDS("output/data.rds")

# Add small floor to w before normalising to avoid log(0) in Stan
data_full$w_norm <- (data_full$w + 1e-6) / mean(data_full$w + 1e-6)

stan_data <- list(
  N     = nrow(data_full),
  logE  = data_full$logE,
  w_norm = data_full$w_norm,
  E_min = 10^min(data_full$logE),
  E_max = 10^max(data_full$logE)
)

rm(data_full)
gc()

# -- Fit ----------------------------------------------------------------------

model <- stan_model("stan/model_spatial.stan")

cat("Starting MCMC: 4 chains x 4500 iter (2 cores)\n")
cat("Start:", format(Sys.time()), "\n\n")

t_start <- Sys.time()

fit_model1 <- sampling(
  model,
  data    = stan_data,
  chains  = 4,
  cores   = 2,
  iter    = 4500,
  warmup  = 1500,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  refresh = 500,
  seed    = 42
)

t_end <- Sys.time()
cat("\n[OK] Done:", format(t_end), "\n")
cat("Elapsed:", round(difftime(t_end, t_start, units = "hours"), 2), "hours\n\n")

saveRDS(fit_model1, "output/fit_model1.rds")
message("Saved: output/fit_model1.rds")

if (SAVE_BACKUP) {
  backup_path <- paste0("output/fit_model1_backup_",
                        format(Sys.time(), "%Y%m%d_%H%M"), ".rds")
  saveRDS(fit_model1, backup_path)
  message("Backup saved: ", backup_path)
}

# -- Diagnostics --------------------------------------------------------------

cat("\n=== CONVERGENCE DIAGNOSTICS ===\n\n")

summary_fit  <- summary(fit_model1)$summary
params       <- c("alpha_astro", "alpha_atmo", "pi_astro")

cat("Parameter estimates:\n")
print(round(summary_fit[params, c("mean", "sd", "n_eff", "Rhat")], 4))
cat("\n")

rhat_vals <- summary_fit[params, "Rhat"]
neff_vals <- summary_fit[params, "n_eff"]

if (all(rhat_vals < 1.01)) {
  cat("[OK] Excellent convergence (Rhat < 1.01)\n")
} else if (all(rhat_vals < 1.05)) {
  cat("[OK] Good convergence (Rhat < 1.05)\n")
} else {
  cat("[WARNING] Rhat > 1.05 for some parameters:\n")
  print(rhat_vals[rhat_vals >= 1.05])
}

cat("Minimum ESS:", min(neff_vals), "\n")

sampler_params <- get_sampler_params(fit_model1, inc_warmup = FALSE)
n_div <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
cat("Divergences:", n_div, "\n")
if (n_div == 0) {
  cat("[OK] No divergences\n\n")
} else if (n_div < 20) {
  cat("[OK] Few divergences — acceptable\n\n")
} else {
  cat("[WARNING] Many divergences — consider raising adapt_delta\n\n")
}

# -- Compute gamma (posterior classification probability) ---------------------

cat("Computing gamma (posterior classification probabilities)...\n")

params_extract <- extract(fit_model1)
alpha_astro    <- mean(params_extract$alpha_astro)
alpha_atmo     <- mean(params_extract$alpha_atmo)
pi_astro       <- mean(params_extract$pi_astro)

# Reload data (was removed from memory above)
data_full        <- readRDS("output/data.rds")
data_full$w_norm <- (data_full$w + 1e-6) / mean(data_full$w + 1e-6)

E     <- 10^data_full$logE
log_E <- log(E)
E_min <- min(E)
E_max <- max(E)

# Truncated power-law normalisations
norm_astro <- (alpha_astro - 1) / (E_min^(1 - alpha_astro) - E_max^(1 - alpha_astro))
norm_atmo  <- (alpha_atmo  - 1) / (E_min^(1 - alpha_atmo)  - E_max^(1 - alpha_atmo))

# Log-likelihoods (spatial weight applies to astrophysical component only)
log_lik_astro <- log(data_full$w_norm) + log(norm_astro) - alpha_astro * log_E
log_lik_atmo  <- log(norm_atmo)                          - alpha_atmo  * log_E

log_prob_astro <- log(pi_astro)       + log_lik_astro
log_prob_atmo  <- log(1 - pi_astro)   + log_lik_atmo

gamma_model1 <- exp(log_prob_astro) / (exp(log_prob_astro) + exp(log_prob_atmo))

saveRDS(gamma_model1, "output/gamma_model1.rds")
message("Saved: output/gamma_model1.rds")

# -- Summary ------------------------------------------------------------------

cat("\n=== RESULTS ===\n\n")
cat(sprintf("Astrophysical fraction (pi_astro): %.4f%%\n", pi_astro * 100))
cat(sprintf("gamma range: [%.6f, %.6f]\n", min(gamma_model1), max(gamma_model1)))
cat(sprintf("High-probability events (top 1%%): %d\n",
            sum(gamma_model1 > quantile(gamma_model1, 0.99))))

cat("\nTop 20 astrophysical candidates:\n")
top_idx <- head(order(gamma_model1, decreasing = TRUE), 20)
for (i in seq_along(top_idx)) {
  idx <- top_idx[i]
  cat(sprintf("%2d. Event %6d: P=%.4f, logE=%.2f\n",
              i, idx, gamma_model1[idx], data_full$logE[idx]))
}

cat("\n[OK] Model 1 fit complete.\n")
