# =============================================================================
# 04_fit_models23.R
# PURPOSE : Fit Model 2 (informed priors) and Model 3 (uniform priors),
#           both energy-only (no spatial weights in the Stan likelihood).
#           Spatial information is combined in a second stage (05_posteriors.R).
#
#           Model 2 priors: alpha_astro ~ N(2.2, 0.5), alpha_atmo ~ N(3.7, 0.3),
#                           pi_astro ~ Beta(2, 35)
#           Model 3 priors: alpha_astro ~ Uniform(1, 4),
#                           alpha_atmo  ~ Uniform(1, 5),
#                           pi_astro    ~ Beta(2, 3)
#
# INPUT   : output/data.rds
#           stan/model_informed.stan
#           stan/model_uniform.stan
# OUTPUT  : output/fit_model2.rds, output/gamma_model2.rds
#           output/fit_model3.rds, output/gamma_model3.rds
# RUNTIME : ~ 7 hours total (two sequential fits, 4 chains each)
# WARNING : This script triggers long MCMC runs. Pre-computed results are
#           already in output/. Run only if you want to reproduce the fits.
# =============================================================================

library(rstan)

options(mc.cores = 4)
rstan_options(auto_write = TRUE)

# -- Prepare data -------------------------------------------------------------

data_full <- readRDS("output/data.rds")

stan_data <- list(
  N     = nrow(data_full),
  logE  = data_full$logE,
  E_min = 10^min(data_full$logE),
  E_max = 10^max(data_full$logE)
)

cat("=== SEQUENTIAL FITS (ENERGY-ONLY, NO SPATIAL WEIGHTS) ===\n\n")

# =============================================================================
# FIT 1: Model 2 — Informed priors
# =============================================================================

cat("--- Model 2: Informed priors ---\n")
cat("Start:", format(Sys.time()), "\n\n")

model2 <- stan_model("stan/model_informed.stan")
t1_start <- Sys.time()

fit_model2 <- sampling(
  model2,
  data    = stan_data,
  chains  = 4,
  cores   = 4,
  iter    = 4000,
  warmup  = 1500,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  refresh = 500,
  seed    = 42
)

t1_end <- Sys.time()
saveRDS(fit_model2, "output/fit_model2.rds")
message("Saved: output/fit_model2.rds")

cat("\n[OK] Model 2 done:", format(t1_end), "\n")
cat("Elapsed:", round(difftime(t1_end, t1_start, units = "hours"), 2), "hours\n\n")

rm(fit_model2, model2)
gc()

# 5-minute pause between fits to allow hardware to cool down
# and memory to be fully released before the second chain run.
cat("Pausing 5 minutes before Model 3 fit...\n\n")
Sys.sleep(300)

# =============================================================================
# FIT 2: Model 3 — Uniform priors
# =============================================================================

cat("--- Model 3: Uniform priors ---\n")
cat("Start:", format(Sys.time()), "\n\n")

model3 <- stan_model("stan/model_uniform.stan")
t2_start <- Sys.time()

fit_model3 <- sampling(
  model3,
  data    = stan_data,
  chains  = 4,
  cores   = 4,
  iter    = 4000,
  warmup  = 1500,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  refresh = 500,
  seed    = 4242
)

t2_end <- Sys.time()
saveRDS(fit_model3, "output/fit_model3.rds")
message("Saved: output/fit_model3.rds")

cat("\n[OK] Model 3 done:", format(t2_end), "\n")
cat("Elapsed:", round(difftime(t2_end, t2_start, units = "hours"), 2), "hours\n")
cat("Total elapsed:", round(difftime(t2_end, t1_start, units = "hours"), 2), "hours\n\n")

# =============================================================================
# DIAGNOSTICS AND GAMMA — Model 2 (Informed priors)
# =============================================================================

cat("\n=== DIAGNOSTICS: MODEL 2 (INFORMED PRIORS) ===\n\n")

fit_model2 <- readRDS("output/fit_model2.rds")
summary_m2 <- summary(fit_model2)$summary
params      <- c("alpha_astro", "alpha_atmo", "pi_astro")

cat("Parameter estimates:\n")
print(round(summary_m2[params, c("mean", "sd", "n_eff", "Rhat")], 4))
cat("\n")

rhat_m2 <- summary_m2[params, "Rhat"]
neff_m2  <- summary_m2[params, "n_eff"]

if (all(rhat_m2 < 1.01)) {
  cat("[OK] Excellent convergence (Rhat < 1.01)\n")
} else if (all(rhat_m2 < 1.05)) {
  cat("[OK] Good convergence (Rhat < 1.05)\n")
} else {
  cat("[WARNING] Rhat > 1.05:\n")
  print(rhat_m2[rhat_m2 >= 1.05])
}
cat("Minimum ESS:", round(min(neff_m2), 0), "\n")

sp_m2  <- get_sampler_params(fit_model2, inc_warmup = FALSE)
n_div2 <- sum(sapply(sp_m2, function(x) sum(x[, "divergent__"])))
cat("Divergences:", n_div2, if (n_div2 == 0) "[OK]\n\n" else
    if (n_div2 < 20) "[OK - few]\n\n" else "[WARNING - many]\n\n")

# Gamma for Model 2 (energy-only, no spatial weight)
cat("Computing gamma for Model 2...\n")
pe2          <- extract(fit_model2)
alpha_astro2 <- mean(pe2$alpha_astro)
alpha_atmo2  <- mean(pe2$alpha_atmo)
pi_astro2    <- mean(pe2$pi_astro)

cat(sprintf("  alpha_astro = %.4f | alpha_atmo = %.4f | pi_astro = %.4f%%\n",
            alpha_astro2, alpha_atmo2, pi_astro2 * 100))

E     <- 10^data_full$logE
log_E <- log(E)
E_min <- min(E); E_max <- max(E)

norm_astro2 <- (alpha_astro2 - 1) / (E_min^(1 - alpha_astro2) - E_max^(1 - alpha_astro2))
norm_atmo2  <- (alpha_atmo2  - 1) / (E_min^(1 - alpha_atmo2)  - E_max^(1 - alpha_atmo2))

log_lk_a2  <- log(norm_astro2) - alpha_astro2 * log_E
log_lk_t2  <- log(norm_atmo2)  - alpha_atmo2  * log_E
gamma_m2   <- exp(log(pi_astro2) + log_lk_a2) /
              (exp(log(pi_astro2) + log_lk_a2) + exp(log(1 - pi_astro2) + log_lk_t2))

saveRDS(gamma_m2, "output/gamma_model2.rds")
message("Saved: output/gamma_model2.rds")

cat("gamma range (M2): [", round(min(gamma_m2), 5), ",",
    round(max(gamma_m2), 5), "]\n")
cat("Top 20 candidates (Model 2):\n")
top2 <- head(order(gamma_m2, decreasing = TRUE, na.last = TRUE), 20)
for (i in seq_along(top2)) {
  idx <- top2[i]
  cat(sprintf("%2d. Event %6d: P=%.4f, logE=%.2f\n",
              i, idx, gamma_m2[idx], data_full$logE[idx]))
}

rm(fit_model2, pe2, sp_m2)
gc()

# =============================================================================
# DIAGNOSTICS AND GAMMA — Model 3 (Uniform priors)
# =============================================================================

cat("\n=== DIAGNOSTICS: MODEL 3 (UNIFORM PRIORS) ===\n\n")

fit_model3 <- readRDS("output/fit_model3.rds")
summary_m3 <- summary(fit_model3)$summary

cat("Parameter estimates:\n")
print(round(summary_m3[params, c("mean", "sd", "n_eff", "Rhat")], 4))
cat("\n")

rhat_m3 <- summary_m3[params, "Rhat"]
neff_m3  <- summary_m3[params, "n_eff"]

if (all(rhat_m3 < 1.01)) {
  cat("[OK] Excellent convergence (Rhat < 1.01)\n")
} else if (all(rhat_m3 < 1.05)) {
  cat("[OK] Good convergence (Rhat < 1.05)\n")
} else {
  cat("[WARNING] Rhat > 1.05:\n")
  print(rhat_m3[rhat_m3 >= 1.05])
}
cat("Minimum ESS:", round(min(neff_m3), 0), "\n")

sp_m3  <- get_sampler_params(fit_model3, inc_warmup = FALSE)
n_div3 <- sum(sapply(sp_m3, function(x) sum(x[, "divergent__"])))
cat("Divergences:", n_div3, if (n_div3 == 0) "[OK]\n\n" else
    if (n_div3 < 20) "[OK - few]\n\n" else "[WARNING - many]\n\n")

# Gamma for Model 3 (energy-only, no spatial weight)
cat("Computing gamma for Model 3...\n")
pe3          <- extract(fit_model3)
alpha_astro3 <- mean(pe3$alpha_astro)
alpha_atmo3  <- mean(pe3$alpha_atmo)
pi_astro3    <- mean(pe3$pi_astro)

cat(sprintf("  alpha_astro = %.4f | alpha_atmo = %.4f | pi_astro = %.4f%%\n",
            alpha_astro3, alpha_atmo3, pi_astro3 * 100))

norm_astro3 <- (alpha_astro3 - 1) / (E_min^(1 - alpha_astro3) - E_max^(1 - alpha_astro3))
norm_atmo3  <- (alpha_atmo3  - 1) / (E_min^(1 - alpha_atmo3)  - E_max^(1 - alpha_atmo3))

log_lk_a3  <- log(norm_astro3) - alpha_astro3 * log_E
log_lk_t3  <- log(norm_atmo3)  - alpha_atmo3  * log_E
gamma_m3   <- exp(log(pi_astro3) + log_lk_a3) /
              (exp(log(pi_astro3) + log_lk_a3) + exp(log(1 - pi_astro3) + log_lk_t3))

saveRDS(gamma_m3, "output/gamma_model3.rds")
message("Saved: output/gamma_model3.rds")

cat("gamma range (M3): [", round(min(gamma_m3), 5), ",",
    round(max(gamma_m3), 5), "]\n")
cat("Top 20 candidates (Model 3):\n")
top3 <- head(order(gamma_m3, decreasing = TRUE, na.last = TRUE), 20)
for (i in seq_along(top3)) {
  idx <- top3[i]
  cat(sprintf("%2d. Event %6d: P=%.4f, logE=%.2f\n",
              i, idx, gamma_m3[idx], data_full$logE[idx]))
}

# =============================================================================
# CROSS-MODEL COMPARISON (Models 2 vs 3)
# =============================================================================

cat("\n=== MODEL 2 vs MODEL 3 PARAMETER COMPARISON ===\n\n")

cat(sprintf("alpha_astro:  %.4f (M2)  vs  %.4f (M3)  [diff: %.4f]\n",
            alpha_astro2, alpha_astro3, abs(alpha_astro2 - alpha_astro3)))
cat(sprintf("alpha_atmo:   %.4f (M2)  vs  %.4f (M3)  [diff: %.4f]\n",
            alpha_atmo2,  alpha_atmo3,  abs(alpha_atmo2  - alpha_atmo3)))
cat(sprintf("pi_astro:     %.4f (M2)  vs  %.4f (M3)  [diff: %.4f]\n\n",
            pi_astro2,    pi_astro3,    abs(pi_astro2    - pi_astro3)))

cat("Gamma correlation (M2 vs M3):",
    round(cor(gamma_m2, gamma_m3, use = "complete.obs"), 6), "\n")
cat("Gamma MAE (M2 vs M3):",
    round(mean(abs(gamma_m2 - gamma_m3), na.rm = TRUE), 6), "\n")

cat("\n[OK] Models 2 and 3 fit complete.\n")
