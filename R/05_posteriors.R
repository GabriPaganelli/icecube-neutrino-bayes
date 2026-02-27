# =============================================================================
# 05_posteriors.R
# PURPOSE : Stage 2 — combine energy-only posteriors (gamma_E) with spatial
#           weights (w) via Bayesian odds updating for Models 2 and 3.
#
#           The update formula is:
#             Odds_E    = gamma_E / (1 - gamma_E)
#             Odds*     = Odds_E * w
#             gamma*    = Odds* / (1 + Odds*)
#
#           Model 1 is excluded here because its spatial weights are already
#           integrated into the Stan likelihood; gamma_model1.rds is final.
#
# INPUT   : output/gamma_model2.rds
#           output/gamma_model3.rds
#           output/data.rds  (provides the w column)
# OUTPUT  : output/gamma_star_model2.rds
#           output/gamma_star_model3.rds
# RUNTIME : < 1 minute
# =============================================================================

# -- Load inputs --------------------------------------------------------------

gamma_E_m2 <- readRDS("output/gamma_model2.rds")
gamma_E_m3 <- readRDS("output/gamma_model3.rds")

data_full  <- readRDS("output/data.rds")
w_space    <- data_full$w

N <- length(gamma_E_m2)
stopifnot(length(gamma_E_m3) == N, length(w_space) == N)

message(sprintf("Loaded %d events. w_space range: [%.4f, %.4f]",
                N, min(w_space), max(w_space)))

# -- Bayesian odds combination ------------------------------------------------

combine_energy_spatial <- function(gamma_E, w_space, epsilon = 1e-10) {
  # Clip gamma_E away from 0 and 1 for numerical stability
  gamma_clipped   <- pmax(pmin(gamma_E, 1 - epsilon), epsilon)
  # Convert to odds, update with spatial weight, convert back
  odds_E          <- gamma_clipped / (1 - gamma_clipped)
  odds_combined   <- odds_E * w_space
  gamma_star      <- odds_combined / (1 + odds_combined)
  return(gamma_star)
}

# -- Model 2 ------------------------------------------------------------------

cat("\n=== MODEL 2: INFORMED PRIORS ===\n")
cat(sprintf("gamma_E range:   [%.6f, %.6f]\n", min(gamma_E_m2), max(gamma_E_m2)))

gamma_star_m2 <- combine_energy_spatial(gamma_E_m2, w_space)

cat(sprintf("gamma* range:    [%.6f, %.6f]\n", min(gamma_star_m2), max(gamma_star_m2)))
cat("gamma* quantiles (M2):\n")
print(round(quantile(gamma_star_m2, probs = c(0.90, 0.95, 0.99, 0.999)), 6))

saveRDS(gamma_star_m2, "output/gamma_star_model2.rds")
message("Saved: output/gamma_star_model2.rds")

# -- Model 3 ------------------------------------------------------------------

cat("\n=== MODEL 3: UNIFORM PRIORS ===\n")
cat(sprintf("gamma_E range:   [%.6f, %.6f]\n", min(gamma_E_m3), max(gamma_E_m3)))

gamma_star_m3 <- combine_energy_spatial(gamma_E_m3, w_space)

cat(sprintf("gamma* range:    [%.6f, %.6f]\n", min(gamma_star_m3), max(gamma_star_m3)))
cat("gamma* quantiles (M3):\n")
print(round(quantile(gamma_star_m3, probs = c(0.90, 0.95, 0.99, 0.999)), 6))

saveRDS(gamma_star_m3, "output/gamma_star_model3.rds")
message("Saved: output/gamma_star_model3.rds")

# -- Impact of spatial information --------------------------------------------

cat("\n=== IMPACT OF SPATIAL INFORMATION ===\n\n")

cat("Model 2 — energy-only vs combined:\n")
cat(sprintf("  Mean gamma_E:  %.6f  |  Mean gamma*: %.6f\n",
            mean(gamma_E_m2), mean(gamma_star_m2)))
cat(sprintf("  Top 1%% cut:   %.6f  ->  %.6f\n",
            quantile(gamma_E_m2, 0.99), quantile(gamma_star_m2, 0.99)))

cat("\nModel 3 — energy-only vs combined:\n")
cat(sprintf("  Mean gamma_E:  %.6f  |  Mean gamma*: %.6f\n",
            mean(gamma_E_m3), mean(gamma_star_m3)))
cat(sprintf("  Top 1%% cut:   %.6f  ->  %.6f\n",
            quantile(gamma_E_m3, 0.99), quantile(gamma_star_m3, 0.99)))

# -- Top 20 candidates --------------------------------------------------------

print_top20 <- function(gamma_star, gamma_E, w, logE, label) {
  cat(sprintf("\n--- Top 20 candidates: %s ---\n", label))
  cat("Rank | Event ID | gamma*   | logE  | gamma_E | w_space\n")
  cat("-----|----------|----------|-------|---------|--------\n")
  top_idx <- head(order(gamma_star, decreasing = TRUE), 20)
  for (i in seq_along(top_idx)) {
    idx <- top_idx[i]
    cat(sprintf("%4d | %8d | %.6f | %5.2f | %.4f  | %.4f\n",
                i, idx, gamma_star[idx], logE[idx], gamma_E[idx], w[idx]))
  }
}

print_top20(gamma_star_m2, gamma_E_m2, w_space, data_full$logE, "Model 2 (Informed)")
print_top20(gamma_star_m3, gamma_E_m3, w_space, data_full$logE, "Model 3 (Uniform)")

cat("\n[OK] Stage 2 posteriors complete.\n")
