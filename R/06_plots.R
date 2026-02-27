# =============================================================================
# 06_plots.R
# PURPOSE : Produce all diagnostic, posterior, and comparison figures.
#           Figures are saved as PDFs to output/figures/.
#
# INPUT   : output/fit_model1.rds, fit_model2.rds, fit_model3.rds
#           output/gamma_model1.rds, gamma_model2.rds, gamma_model3.rds
#           output/gamma_star_model2.rds, gamma_star_model3.rds
#           output/data.rds, output/sources.rds
#
# OUTPUT  : output/figures/map.pdf               — sky map (Model 1)
#           output/figures/diagnostics.pdf        — MCMC trace/density/ACF
#           output/figures/posteriors.pdf         — parameter comparison
#           output/figures/ppc_model1.pdf         — posterior predictive checks
#           output/figures/gamma_comparison.pdf   — final gamma* across models
#
# RUNTIME : < 5 minutes
# =============================================================================

library(rstan)
library(bayesplot)
library(ggplot2)
library(gridExtra)
library(posterior)
library(patchwork)
library(ggnewscale)

if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)

# -- Load all results ---------------------------------------------------------

data    <- readRDS("output/data.rds")
sources <- readRDS("output/sources.rds")

fit1 <- readRDS("output/fit_model1.rds")
fit2 <- readRDS("output/fit_model2.rds")
fit3 <- readRDS("output/fit_model3.rds")

gamma_m1      <- readRDS("output/gamma_model1.rds")
gamma_m2      <- readRDS("output/gamma_model2.rds")
gamma_m3      <- readRDS("output/gamma_model3.rds")
gamma_star_m2 <- readRDS("output/gamma_star_model2.rds")
gamma_star_m3 <- readRDS("output/gamma_star_model3.rds")

# =============================================================================
# FIGURE 1: Sky map — high-probability events overlaid on source positions
# =============================================================================

message("Generating sky map...")

threshold    <- quantile(gamma_m1, 0.90)
high_idx     <- which(gamma_m1 > threshold)
data_high    <- data[high_idx, ]
gamma_high   <- gamma_m1[high_idx]

p_map <- ggplot() +
  geom_point(
    data = data.frame(ra = data_high$ra, dec = data_high$dec, gamma = gamma_high),
    aes(x = ra, y = dec, color = gamma, size = gamma),
    alpha = 0.6, shape = 16) +
  scale_color_gradient(low = "lightblue", high = "darkblue",
                       name = expression(gamma)) +
  scale_size_continuous(range = c(0.5, 3), guide = "none") +
  ggnewscale::new_scale_color() +
  geom_point(
    data = sources,
    aes(x = RA_deg, y = Dec_deg, fill = ns_hat),
    color = "black", shape = 21, size = 4, stroke = 1.5, alpha = 0.8) +
  scale_fill_gradient(low = "yellow", high = "red",
                      name = expression(hat(n)[s])) +
  labs(x = "RA (deg)", y = "Dec (deg)",
       title = "Sources and high-probability astrophysical events (Model 1)",
       subtitle = paste0("Top 10%  (gamma > ", round(threshold, 3), ")")) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        panel.grid.minor = element_line(linetype = "dotted", color = "gray90"))

pdf("output/figures/map.pdf", width = 10, height = 5)
print(p_map)
dev.off()
message("Saved: output/figures/map.pdf")

# =============================================================================
# FIGURE 2: MCMC diagnostics — trace, density, ACF for all three models
# =============================================================================

message("Generating diagnostics...")

diagnostic_plots <- function(fit, model_name) {
  p_trace <- mcmc_trace(fit, pars = c("alpha_astro", "alpha_atmo", "pi_astro")) +
    ggtitle(paste(model_name, "— Trace Plots"))
  p_dens  <- mcmc_dens_overlay(fit, pars = c("alpha_astro", "alpha_atmo", "pi_astro")) +
    ggtitle(paste(model_name, "— Posterior Densities"))
  p_acf   <- mcmc_acf(fit, pars = c("alpha_astro", "alpha_atmo", "pi_astro")) +
    ggtitle(paste(model_name, "— Autocorrelation"))
  list(trace = p_trace, density = p_dens, acf = p_acf)
}

check_diagnostics <- function(fit, model_name) {
  summ <- summary(fit)$summary
  cat("\n===", model_name, "===\n")
  cat("ESS (n_eff):\n"); print(round(summ[, "n_eff"], 0))
  cat("Rhat:\n"); print(round(summ[, "Rhat"], 4))
  if (any(summ[, "Rhat"] > 1.01, na.rm = TRUE))
    warning(paste(model_name, ": some Rhat > 1.01"))
  if (any(summ[, "n_eff"] < 400, na.rm = TRUE))
    warning(paste(model_name, ": some ESS < 400"))
}

check_diagnostics(fit1, "Model 1 (Spatial)")
check_diagnostics(fit2, "Model 2 (Informed)")
check_diagnostics(fit3, "Model 3 (Uniform)")

diag1 <- diagnostic_plots(fit1, "Model 1 (Spatial)")
diag2 <- diagnostic_plots(fit2, "Model 2 (Informed)")
diag3 <- diagnostic_plots(fit3, "Model 3 (Uniform)")

pdf("output/figures/diagnostics.pdf", width = 12, height = 8)
for (d in list(diag1, diag2, diag3)) {
  print(grid.arrange(d$trace, d$density, d$acf, ncol = 1))
}
dev.off()
message("Saved: output/figures/diagnostics.pdf")

# =============================================================================
# FIGURE 3: Posterior parameter comparison — all three models
# =============================================================================

message("Generating parameter comparison plots...")

# Prior-posterior overlay for pi_astro (Model 1)
samples_pi1 <- as.vector(as.matrix(fit1, pars = "pi_astro"))
x_prior     <- seq(0, 0.1, length.out = 1000)
df_prior    <- data.frame(x = x_prior, density = dbeta(x_prior, 1, 99))

p_prior_post <- ggplot() +
  geom_density(data = data.frame(pi_astro = samples_pi1),
               aes(x = pi_astro),
               fill = "steelblue", alpha = 0.3, linewidth = 1) +
  geom_line(data = df_prior, aes(x = x, y = density),
            linewidth = 1, linetype = "dashed") +
  labs(title = "Posterior vs Prior for pi_astro (Model 1)",
       y = "Density", x = expression(pi[astro])) +
  coord_cartesian(xlim = c(0, 0.1)) +
  theme_minimal()

# Model comparison: density overlay per parameter
compare_models_facet <- function(fit1, fit2, fit3) {
  params <- c("alpha_astro", "alpha_atmo", "pi_astro")
  df <- do.call(rbind, lapply(params, function(par) {
    s1 <- as.vector(as.matrix(fit1, pars = par))
    s2 <- as.vector(as.matrix(fit2, pars = par))
    s3 <- as.vector(as.matrix(fit3, pars = par))
    data.frame(
      value     = c(s1, s2, s3),
      model     = rep(c("M1 (Spatial)", "M2 (Informed)", "M3 (Uniform)"),
                      c(length(s1), length(s2), length(s3))),
      parameter = par
    )
  }))
  ggplot(df, aes(x = value, color = model, fill = model)) +
    geom_density(alpha = 0.3, linewidth = 1) +
    facet_wrap(~ parameter, scales = "free") +
    ggtitle("Posterior Comparison — All Parameters") +
    theme_minimal() +
    theme(legend.position = "bottom")
}

compare_summary <- function(fit1, fit2, fit3, par) {
  s1 <- as.vector(as.matrix(fit1, pars = par))
  s2 <- as.vector(as.matrix(fit2, pars = par))
  s3 <- as.vector(as.matrix(fit3, pars = par))
  cat(sprintf("\n[%s]\n", par))
  for (nm in c("M1", "M2", "M3")) {
    s <- switch(nm, M1 = s1, M2 = s2, M3 = s3)
    cat(sprintf("  %s: mean=%.4f  sd=%.4f  [%.4f, %.4f]\n",
                nm, mean(s), sd(s),
                quantile(s, 0.025), quantile(s, 0.975)))
  }
}

p_compare <- compare_models_facet(fit1, fit2, fit3)

cat("\n=== PARAMETER SUMMARY ===\n")
for (par in c("alpha_astro", "alpha_atmo", "pi_astro"))
  compare_summary(fit1, fit2, fit3, par)

cat("\n=== FULL POSTERIOR SUMMARIES ===\n")
cat("\n-- Model 1 --\n"); print(summary(fit1)$summary)
cat("\n-- Model 2 --\n"); print(summary(fit2)$summary)
cat("\n-- Model 3 --\n"); print(summary(fit3)$summary)

pdf("output/figures/posteriors.pdf", width = 12, height = 8)
print(p_prior_post)
print(p_compare)
dev.off()
message("Saved: output/figures/posteriors.pdf")

# =============================================================================
# FIGURE 4: Posterior predictive checks — Model 1
# =============================================================================

message("Generating posterior predictive checks (Model 1)...")

generate_truncated_powerlaw <- function(alpha, E_min, E_max, n) {
  u <- runif(n)
  if (abs(alpha - 1) < 1e-6) {
    E_min * exp(u * log(E_max / E_min))
  } else {
    E_min_pow <- E_min^(1 - alpha)
    E_max_pow <- E_max^(1 - alpha)
    (E_min_pow + u * (E_max_pow - E_min_pow))^(1 / (1 - alpha))
  }
}

posterior_samples <- extract(fit1)
E_min    <- 10^min(data$logE)
E_max    <- 10^max(data$logE)
n_events <- length(data$logE)

# Generate 200 full predictive draws (subsampled for heavy plots)
n_sims_full  <- 200
post_indices <- sample(length(posterior_samples$alpha_astro), n_sims_full)
y_rep_full   <- matrix(NA, nrow = n_sims_full, ncol = n_events)

cat("Generating", n_sims_full, "posterior predictive samples...\n")
for (i in seq_len(n_sims_full)) {
  idx       <- post_indices[i]
  a_astro   <- posterior_samples$alpha_astro[idx]
  a_atmo    <- posterior_samples$alpha_atmo[idx]
  pi_a      <- posterior_samples$pi_astro[idx]
  n_astro   <- rbinom(1, n_events, pi_a)
  n_atmo    <- n_events - n_astro
  E_astro   <- generate_truncated_powerlaw(a_astro, E_min, E_max, n_astro)
  E_atmo    <- generate_truncated_powerlaw(a_atmo,  E_min, E_max, n_atmo)
  y_rep_full[i, ] <- log10(c(E_astro, E_atmo))
  if (i %% 50 == 0) { cat("  ", i, "/", n_sims_full, "\n"); gc() }
}

# Subsample 10k events for lighter plots
n_sub       <- 10000
sub_idx     <- sample(n_events, n_sub)
y_rep_sub   <- y_rep_full[, sub_idx]
data_sub    <- data$logE[sub_idx]

p_ppc1 <- ppc_dens_overlay(data_sub, y_rep_sub[1:30, ]) +
  labs(title = "PPC: Energy Distribution (Model 1)",
       x = expression(log[10](E / GeV)), y = "Density") +
  theme_minimal()

p_ppc2 <- ppc_hist(data_sub, y_rep_sub[1:8, ], binwidth = 0.15) +
  labs(title = "PPC: Histogram Grid (Model 1)") +
  theme_minimal()

p_ppc3 <- grid.arrange(
  ppc_stat(data$logE, y_rep_full, stat = "mean") +
    labs(title = "PPC: Mean Energy") + theme_minimal(),
  ppc_stat(data$logE, y_rep_full, stat = "sd") +
    labs(title = "PPC: SD Energy") + theme_minimal(),
  ppc_stat(data$logE, y_rep_full, stat = "median") +
    labs(title = "PPC: Median Energy") + theme_minimal(),
  ncol = 3
)

# Aggregated ribbon (90% predictive interval)
bins           <- seq(min(data$logE), max(data$logE), length.out = 100)
bin_centers    <- (bins[-1] + bins[-length(bins)]) / 2
density_matrix <- matrix(NA, nrow = n_sims_full, ncol = length(bin_centers))
for (i in seq_len(n_sims_full))
  density_matrix[i, ] <- hist(y_rep_full[i, ], breaks = bins, plot = FALSE)$density

plot_df <- data.frame(
  x      = bin_centers,
  median = apply(density_matrix, 2, median),
  lower  = apply(density_matrix, 2, quantile, 0.05),
  upper  = apply(density_matrix, 2, quantile, 0.95)
)

p_ppc4 <- ggplot() +
  geom_ribbon(data = plot_df, aes(x = x, ymin = lower, ymax = upper),
              fill = "#3182bd", alpha = 0.3) +
  geom_line(data = plot_df, aes(x = x, y = median),
            color = "#3182bd", linewidth = 1) +
  geom_density(data = data.frame(E = data$logE), aes(x = E),
               color = "black", linewidth = 1.2) +
  labs(title = "Aggregated PPC: 90% Predictive Interval (Model 1)",
       subtitle = "Blue: predictive interval (200 draws)  |  Black: observed",
       x = expression(log[10](E / GeV)), y = "Density") +
  theme_minimal()

pdf("output/figures/ppc_model1.pdf", width = 12, height = 8)
print(p_ppc1)
print(p_ppc2)
print(p_ppc3)
print(p_ppc4)
dev.off()
message("Saved: output/figures/ppc_model1.pdf")

rm(y_rep_full, y_rep_sub, density_matrix)
gc()

# =============================================================================
# FIGURE 5: Final gamma* comparison across all three models
# =============================================================================
# gamma_m1        = Model 1 final posterior (spatial integrated in Stan)
# gamma_star_m2   = Model 2 final posterior (energy-only + post-hoc spatial)
# gamma_star_m3   = Model 3 final posterior (energy-only + post-hoc spatial)
# Note: max observed gamma ~ 0.37; threshold-based counts are not meaningful.
# =============================================================================

message("Generating gamma* comparison plots...")

df_gamma <- rbind(
  data.frame(gamma = gamma_m1,      model = "M1 (Spatial)",   type = "final"),
  data.frame(gamma = gamma_star_m2, model = "M2 (Informed*)", type = "final"),
  data.frame(gamma = gamma_star_m3, model = "M3 (Uniform*)",  type = "final")
)

# Density overlay of final gamma* distributions
p_g1 <- ggplot(df_gamma, aes(x = gamma, color = model, fill = model)) +
  geom_density(alpha = 0.25, linewidth = 1) +
  labs(title = "Final classification probability: all models",
       subtitle = "M1: spatial in Stan | M2/M3: post-hoc spatial combination",
       x = expression(gamma^"*"), y = "Density") +
  theme_minimal() +
  theme(legend.position = "bottom")

# High-end density (zoom in to top 1%)
cut_01pct <- quantile(c(gamma_m1, gamma_star_m2, gamma_star_m3), 0.99)
p_g2 <- ggplot(subset(df_gamma, gamma > cut_01pct),
               aes(x = gamma, color = model, fill = model)) +
  geom_density(alpha = 0.25, linewidth = 1) +
  labs(title = "Top 1% classification probability — distribution tail",
       x = expression(gamma^"*"), y = "Density") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Quantile summary table
cat("\n=== GAMMA* QUANTILE SUMMARY ===\n")
probs <- c(0.90, 0.95, 0.99, 0.999)
for (nm in c("M1", "M2*", "M3*")) {
  g <- switch(nm, "M1" = gamma_m1, "M2*" = gamma_star_m2, "M3*" = gamma_star_m3)
  cat(sprintf("\n%s: %s\n", nm,
              paste(sprintf("%.6f", quantile(g, probs)),
                    paste0("(", probs * 100, "%)"), collapse = " | ")))
}

# Top 20 candidates — all three models
print_top20_comparison <- function(gamma_vec, logE_vec, label) {
  top <- head(order(gamma_vec, decreasing = TRUE), 20)
  cat(sprintf("\n--- Top 20: %s ---\n", label))
  cat("Rank | Event ID | gamma*   | logE\n")
  cat("-----|----------|----------|-----\n")
  for (i in seq_along(top)) {
    idx <- top[i]
    cat(sprintf("%4d | %8d | %.6f | %.2f\n",
                i, idx, gamma_vec[idx], logE_vec[idx]))
  }
}

print_top20_comparison(gamma_m1,      data$logE, "Model 1 (Spatial)")
print_top20_comparison(gamma_star_m2, data$logE, "Model 2 (Informed + Spatial*)")
print_top20_comparison(gamma_star_m3, data$logE, "Model 3 (Uniform + Spatial*)")

# Cross-model gamma* scatter (M1 vs M2* and M1 vs M3*)
df_scatter_sub <- data.frame(
  gamma_m1      = gamma_m1[seq(1, length(gamma_m1), by = 100)],
  gamma_star_m2 = gamma_star_m2[seq(1, length(gamma_m1), by = 100)],
  gamma_star_m3 = gamma_star_m3[seq(1, length(gamma_m1), by = 100)]
)

p_g3 <- ggplot(df_scatter_sub, aes(x = gamma_m1, y = gamma_star_m2)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "gamma* comparison: M1 vs M2*",
       x = "M1 (Spatial, integrated)", y = "M2 (Informed + post-hoc spatial)") +
  theme_minimal()

p_g4 <- ggplot(df_scatter_sub, aes(x = gamma_m1, y = gamma_star_m3)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "gamma* comparison: M1 vs M3*",
       x = "M1 (Spatial, integrated)", y = "M3 (Uniform + post-hoc spatial)") +
  theme_minimal()

pdf("output/figures/gamma_comparison.pdf", width = 12, height = 8)
print(p_g1)
print(p_g2)
print(grid.arrange(p_g3, p_g4, ncol = 2))
dev.off()
message("Saved: output/figures/gamma_comparison.pdf")

message("\n[OK] All figures saved to output/figures/")
