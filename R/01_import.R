# =============================================================================
# 01_import.R
# PURPOSE : Load raw IceCube CSV event files and the source catalogue.
#           Combine all seasons into a single dataset.
# INPUT   : data/raw/IC86_*.csv, data/sources.txt
# OUTPUT  : output/icecube_all.rds
# RUNTIME : < 1 minute
# =============================================================================

library(tidyverse)

# -- Load event files ---------------------------------------------------------

data_dir <- "data/raw"
files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)

read_icecube <- function(f) {
  cat("Loading:", basename(f), "\n")
  df <- read.table(
    f,
    header = TRUE,
    comment.char = "#",
    stringsAsFactors = FALSE
  )
  colnames(df) <- c("mjd", "logE", "ang_err", "ra", "dec", "azimuth", "zenith")
  df$season <- gsub("_exp.*", "", basename(f))
  return(df)
}

data_list <- lapply(files, read_icecube)
icecube_all <- do.call(rbind, data_list)

cat(sprintf("\n[OK] Loaded %d events from %d seasons\n",
            nrow(icecube_all), length(files)))

# -- Load source catalogue ----------------------------------------------------

sources_raw <- read.table("data/sources.txt", sep = ",", header = TRUE)
cat(sprintf("[OK] Source catalogue: %d sources\n", nrow(sources_raw)))

# -- Save ---------------------------------------------------------------------

saveRDS(icecube_all, "output/icecube_all.rds")
message("Saved: output/icecube_all.rds")
