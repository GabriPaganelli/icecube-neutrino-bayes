# =============================================================================
# 02_preprocess.R
# PURPOSE : Clean the raw event dataset, run exploratory summaries,
#           compute Gaussian spatial weights for each event relative to
#           110 known point sources.
# INPUT   : output/icecube_all.rds, data/sources.txt
# OUTPUT  : output/data.rds     — cleaned events + spatial weight column
#           output/sources.rds  — source catalogue with sigma_src2 column
# RUNTIME : ~ 5 minutes (spatial weight matrix: N × J = 880k × 110)
# NOTE    : theta_mat and w_mat are large intermediate matrices (~1 GB total)
#           and are intentionally NOT saved to disk; only data$w is retained.
# =============================================================================

library(tidyverse)
library(patchwork)

# =============================================================================
# 1. LOAD DATA
# =============================================================================

data    <- readRDS("output/icecube_all.rds")
sources <- read.table("data/sources.txt", sep = ",", header = TRUE)

cat(sprintf("Events loaded: %d\n", nrow(data)))
cat(sprintf("Sources loaded: %d\n", nrow(sources)))

# =============================================================================
# 2. EDA UTILITY FUNCTIONS
# (exploratory; safe to skip — none of these produce pipeline outputs)
# =============================================================================

# Summarise each column: type, range, NA count, unique count
summarize_variables <- function(data) {

  if (!is.data.frame(data) && !inherits(data, c("tbl_df", "data.table"))) {
    stop("Input must be a data.frame, tibble, or data.table")
  }

  if (ncol(data) == 0) {
    return(data.frame(
      variable = character(0), type = character(0),
      range    = character(0), n_na = integer(0), n_unique = integer(0),
      stringsAsFactors = FALSE
    ))
  }

  n_cols <- ncol(data)
  result <- data.frame(
    variable = character(n_cols), type     = character(n_cols),
    range    = character(n_cols), n_na     = integer(n_cols),
    n_unique = integer(n_cols),   stringsAsFactors = FALSE
  )

  compute_range_string <- function(col_data, n_unique, n_na) {
    if (n_na == length(col_data)) return("all NA")
    if (is.numeric(col_data)) {
      return(paste0("[", min(col_data, na.rm = TRUE), " - ",
                        max(col_data, na.rm = TRUE), "]"))
    } else if (is.factor(col_data) || is.character(col_data)) {
      return(paste0(n_unique, " categories"))
    } else if (is.logical(col_data)) {
      return(paste(unique(col_data[!is.na(col_data)]), collapse = ", "))
    } else if (inherits(col_data, "Date")) {
      mn <- min(col_data, na.rm = TRUE); mx <- max(col_data, na.rm = TRUE)
      if (is.na(mn) || is.na(mx)) return("all NA")
      return(paste0("[", mn, " - ", mx, "]"))
    } else {
      return(paste0(n_unique, " unique values"))
    }
  }

  for (i in seq_len(ncol(data))) {
    col_data <- data[[i]]
    n_na     <- sum(is.na(col_data))
    n_unique <- length(unique(col_data[!is.na(col_data)]))
    result[i, ] <- list(
      variable = names(data)[i],
      type     = class(col_data)[1],
      range    = compute_range_string(col_data, n_unique, n_na),
      n_na     = n_na,
      n_unique = n_unique
    )
  }
  return(result)
}

# Correlation matrix with fallback rendering
plot_corrplot <- function(cor_matrix, title, method) {
  if (requireNamespace("corrplot", quietly = TRUE)) {
    tryCatch({
      corrplot::corrplot(cor_matrix,
                         method   = "color", type  = "upper",
                         order    = "original", tl.col = "black",
                         tl.srt   = 45, diag = FALSE,
                         title    = paste(title, "-", method),
                         mar      = c(0, 0, 2, 0),
                         addCoef.col = if (ncol(cor_matrix) <= 10) "black" else NULL)
    }, error = function(e) {
      warning("corrplot failed, falling back to image plot")
      plot_imageplot_matrix(cor_matrix, title, method)
    })
  } else {
    warning("Package corrplot not available, using image plot")
    plot_imageplot_matrix(cor_matrix, title, method)
  }
}

plot_imageplot_matrix <- function(cor_matrix, title, method) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  p          <- ncol(cor_matrix)
  var_names  <- colnames(cor_matrix)
  cex_size   <- max(0.3, min(1, 15 / p))
  bot_margin <- max(3, min(15, max(nchar(var_names)) * 0.3))
  lft_margin <- max(3, min(15, max(nchar(var_names)) * 0.3))
  par(mar = c(bot_margin, lft_margin, 4, 4))

  image(1:p, 1:p, cor_matrix,
        col   = colorRampPalette(c("red", "white", "blue"))(100),
        zlim  = c(-1, 1), xlab = "", ylab = "",
        main  = paste(title, "\n(Method:", method, ")"),
        axes  = FALSE)

  if (p <= 50) {
    axis(1, at = 1:p, labels = var_names, las = 2, cex.axis = cex_size)
    axis(2, at = 1:p, labels = var_names, las = 1, cex.axis = cex_size)
  } else {
    tick_pos    <- seq(1, p, length.out = min(20, p))
    tick_labels <- var_names[round(tick_pos)]
    axis(1, at = tick_pos, labels = tick_labels, las = 2, cex.axis = 0.3)
    axis(2, at = tick_pos, labels = tick_labels, las = 1, cex.axis = 0.3)
  }

  if (p <= 20) {
    abline(h = 1:p + 0.5, col = "gray", lwd = 0.5)
    abline(v = 1:p + 0.5, col = "gray", lwd = 0.5)
  }
}

# Compute and optionally plot the correlation matrix of numeric columns
plot_correlation_matrix <- function(data, y_col = NULL, method = "pearson",
                                    plot = TRUE,
                                    title = "Correlation Matrix",
                                    p_threshold = 20) {

  if (!is.data.frame(data) && !inherits(data, c("tbl_df", "data.table"))) {
    stop("Input must be a data.frame, tibble, or data.table")
  }
  if (!is.null(y_col) && !y_col %in% names(data)) {
    warning(paste("Column", y_col, "not found; ignoring y_col"))
    y_col <- NULL
  }

  numeric_vars <- sapply(data, is.numeric)
  x_vars <- if (!is.null(y_col)) {
    names(data)[numeric_vars & names(data) != y_col]
  } else {
    names(data)[numeric_vars]
  }

  if (length(x_vars) < 2)
    stop("Need at least 2 numeric columns to compute a correlation matrix")

  x_data    <- data[, x_vars, drop = FALSE]
  n_complete <- sum(complete.cases(x_data))
  if (n_complete < 2)
    stop("Not enough complete observations to compute correlations")

  cor_matrix <- cor(x_data, method = method, use = "pairwise.complete.obs")

  if (any(is.na(diag(cor_matrix))))
    warning("Some variables have zero variance or computation issues")

  if (plot && !is.null(cor_matrix)) {
    p <- ncol(cor_matrix)
    if (p <= p_threshold) {
      plot_corrplot(cor_matrix, title, method)
    } else {
      plot_imageplot_matrix(cor_matrix, title, method)
    }
  }
  return(cor_matrix)
}

# Plot distribution for each column (histogram / bar chart)
plot_variable_distributions <- function(df) {

  if (!is.data.frame(df) && !inherits(df, c("tbl_df", "tbl", "data.table")))
    stop("Input must be a data.frame, tibble, or data.table")

  df <- as.data.frame(df)
  if (nrow(df) == 0 || ncol(df) == 0) stop("data.frame is empty")

  if (!requireNamespace("ggplot2",   quietly = TRUE)) stop("Install ggplot2")
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("Install patchwork")

  plot_list <- list()

  for (col_name in names(df)) {
    col_data <- df[[col_name]]
    if (all(is.na(col_data))) {
      message(paste("Skipping", col_name, ": all values are NA"))
      next
    }

    p <- NULL

    if (is.numeric(col_data)) {
      non_na_data <- col_data[!is.na(col_data)]
      if (length(non_na_data) == 0) next
      n_bins <- min(30, max(10, length(non_na_data) / 10))
      p <- ggplot(df, aes(x = .data[[col_name]])) +
        geom_histogram(fill = "skyblue", color = "black",
                       bins = n_bins, alpha = 0.7) +
        labs(title = paste("Distribution:", col_name), x = col_name,
             y = "Count",
             subtitle = if (sum(is.na(col_data)) > 0)
               paste("NA values:", sum(is.na(col_data))) else NULL) +
        theme_minimal()

    } else if (is.factor(col_data) || is.character(col_data) || is.logical(col_data)) {
      tmp <- df
      tmp[[col_name]] <- as.character(tmp[[col_name]])
      tmp[[col_name]][is.na(tmp[[col_name]])] <- "NAs"
      orig_levels <- if (is.factor(col_data)) levels(col_data) else unique(na.omit(col_data))
      all_levels  <- c(orig_levels, "NAs")
      tmp[[col_name]] <- factor(tmp[[col_name]], levels = all_levels)

      max_cats <- 20
      if (length(levels(tmp[[col_name]])) > max_cats) {
        freq   <- table(tmp[[col_name]])
        top_cats <- names(sort(freq[names(freq) != "NAs"], decreasing = TRUE))[1:(max_cats - 2)]
        tmp[[col_name]] <- as.character(tmp[[col_name]])
        tmp[[col_name]][!tmp[[col_name]] %in% c(top_cats, "NAs")] <- "Other"
        final_levels <- c(top_cats, "Other",
                          if (sum(tmp[[col_name]] == "NAs") > 0) "NAs" else NULL)
        tmp[[col_name]] <- factor(tmp[[col_name]], levels = final_levels)
      }

      p <- ggplot(tmp, aes(x = .data[[col_name]])) +
        geom_bar(fill = "lightgreen", color = "black", alpha = 0.7) +
        labs(title = paste("Distribution:", col_name),
             x = col_name, y = "Count") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

    } else {
      message(paste("Skipping", col_name, ": unsupported type", class(col_data)[1]))
      next
    }

    if (!is.null(p)) plot_list[[col_name]] <- p
  }

  if (length(plot_list) == 0) {
    message("No plots generated.")
    return(invisible(NULL))
  }

  n_plots <- length(plot_list)
  if (n_plots <= 4) {
    result <- wrap_plots(plot_list, ncol = 2)
  } else if (n_plots <= 9) {
    result <- wrap_plots(plot_list, ncol = 3)
  } else {
    plots_per_page <- 9
    n_pages <- ceiling(n_plots / plots_per_page)
    message(paste("Generating", n_pages, "pages of plots"))
    pages <- lapply(seq_len(n_pages), function(pg) {
      idx <- ((pg - 1) * plots_per_page + 1):min(pg * plots_per_page, n_plots)
      wrap_plots(plot_list[idx], ncol = 3)
    })
    result <- pages[[1]]
    if (n_pages > 1) for (i in 2:n_pages) print(pages[[i]])
  }

  return(result)
}

# =============================================================================
# 3. EDA  (exploratory — safe to skip when re-running the pipeline)
# =============================================================================

cat("\n-- Variable summary (full dataset) --\n")
print(summarize_variables(data))
cat("\n-- Variable summary (sources) --\n")
print(summarize_variables(sources))

plot_correlation_matrix(data, title = "Event variable correlations")
print(plot_variable_distributions(data))

# Angular error distribution
cat("\nAngular error (ang_err) quantiles:\n")
print(summary(data$ang_err))
cat("Fraction > 2 deg:", mean(data$ang_err > 2), "\n")
cat("Fraction > 4 deg:", mean(data$ang_err > 4), "\n")

cat("\nSource catalogue class breakdown:\n")
print(table(sources$Class))
plot(sources$RA_deg, sources$Dec_deg,
     xlab = "RA (deg)", ylab = "Dec (deg)",
     main = "Source positions", pch = 19, cex = 0.7)

# =============================================================================
# 4. FILTERING
# =============================================================================

# Remove events with large angular uncertainty (> 4 deg) and local coordinates
data_clean <- data[data$ang_err < 4, -c(6, 7)]   # drop azimuth, zenith
cat(sprintf("\n[OK] After filtering: %d events (removed %d with ang_err >= 4)\n",
            nrow(data_clean), nrow(data) - nrow(data_clean)))

# =============================================================================
# 5. FEATURE ENGINEERING: SPATIAL WEIGHTS
# =============================================================================
# For each event i and source j, compute the Gaussian spatial weight:
#
#   w_{ij} = exp( -theta_{ij}^2 / (2 * sigma_eff^2) )
#
# where sigma_eff^2 = sigma_meas^2 + sigma_src^2,
#   sigma_meas = ang_err  (event angular uncertainty, degrees)
#   sigma_src  = a / (ns_hat + b)^sigma_power  (source size proxy)
#
# The per-event weight is w_i = sum_j w_{ij}, normalised to [0, 1].

# Hyperparameters (spatial weight kernel)
a           <- 0.5    # angular scale (deg)
b           <- 1      # stabiliser so ns_hat = 0 doesn't blow up
sigma_power <- 0.5    # exponent in the source-size power law

sigma_min <- 0.05     # floor for sigma_src (deg)

# Source uncertainty proxy
sources$sigma_src2 <- pmax(sigma_min, a / (sources$ns_hat + b)^sigma_power)^2

# Event measurement variance
data_clean$sigma_meas2 <- data_clean$ang_err^2

# Angular distance matrix (N × J) using the haversine-like formula
deg2rad <- pi / 180
ra_e    <- data_clean$ra      * deg2rad
dec_e   <- data_clean$dec     * deg2rad
ra_s    <- sources$RA_deg     * deg2rad
dec_s   <- sources$Dec_deg    * deg2rad

N <- nrow(data_clean)
J <- nrow(sources)

cat(sprintf("\nComputing angular distance matrix (%d × %d)...\n", N, J))
theta_mat <- matrix(0.0, nrow = N, ncol = J)
for (j in seq_len(J)) {
  cos_d <- sin(dec_e) * sin(dec_s[j]) +
            cos(dec_e) * cos(dec_s[j]) * cos(abs(ra_e - ra_s[j]))
  cos_d <- pmin(1, pmax(-1, cos_d))          # numerical clamp
  theta_mat[, j] <- acos(cos_d) / deg2rad    # degrees
}

cat("Computing spatial weight matrix...\n")
w_mat <- matrix(0.0, nrow = N, ncol = J)
for (j in seq_len(J)) {
  sigma_eff2     <- data_clean$sigma_meas2 + sources$sigma_src2[j]
  w_mat[, j]     <- exp(-theta_mat[, j]^2 / (2 * sigma_eff2))
}

# Aggregate: per-event weight = sum of source contributions, normalised
data_clean$w <- rowSums(w_mat)
data_clean$w <- data_clean$w / max(data_clean$w)

# theta_mat and w_mat are large (~800 MB total) and not needed downstream;
# let them be garbage-collected rather than written to disk.
rm(theta_mat, w_mat)
gc()

cat("[OK] Spatial weights computed. Range: [",
    round(min(data_clean$w), 5), ",", round(max(data_clean$w), 5), "]\n")

# =============================================================================
# 6. SAVE
# =============================================================================

saveRDS(data_clean, "output/data.rds")
message("Saved: output/data.rds")

saveRDS(sources, "output/sources.rds")
message("Saved: output/sources.rds")
