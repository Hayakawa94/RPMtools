# Install and load necessary packages
required_packages <- c(
  'tidyverse', 'odbc', 'dbplyr', 'data.table', 'CatEncoders', 'glue', 'plotly',
  'htmltools', 'dplyr', 'sf', 'gridExtra', 'tidyr', 'lubridate', 'reshape2',
  'reticulate', 'ggplot2', 'ParBayesianOptimization', 'mlbench', 'resample',
  'xgboost', 'Matrix', 'pracma', 'RColorBrewer', 'cartogram', 'tmap', 'spdep',
  'deldir', 'sp', 'purrr', 'DescTools', 'readxl', 'openxlsx', 'fastglm',
  'dtplyr', 'pbapply', 'patchwork', 'shiny', 'writexl'
)

new_packages <- setdiff(required_packages, installed.packages()[, "Package"])
if (length(new_packages) > 0) install.packages(new_packages)

invisible(lapply(required_packages, library, character.only = TRUE))

# Options and utility functions
options(scipen = 999)

KT_create_equal_bin <- function(weight, nbin) {
  cumulative_sum <- cumsum(weight)
  bins <- cut(cumulative_sum, breaks = nbin, labels = FALSE)
  return(bins)
}

# Gini calculations
KT_calc_gini <- function(actual, weight, predicted) {
  df <- data.frame(actual = as.numeric(actual), weight = as.numeric(weight), predicted)
  sorted_idx <- order(df$predicted)
  w_s <- df$weight[sorted_idx]
  a_s <- df$actual[sorted_idx]
  a_c <- cumsum(a_s * w_s)
  w_c <- cumsum(w_s)
  return(1 - 2 * pracma::trapz(w_c / max(w_c), a_c / max(a_c)))
}

KT_calc_gini_norm <- function(actual, weight, predicted) {
  return(KT_calc_gini(actual, weight, predicted) / KT_calc_gini(actual, weight, actual))
}

# Resampling functions
KT_resample_gini <- function(n, actual, weight, predicted, normalize = FALSE) {
  gini_vector <- numeric(n)
  df <- data.frame(actual, weight, predicted)
  set.seed(123)
  for (i in seq_len(n)) {
    sampled_df <- df[sample(nrow(df), replace = TRUE), ]
    gini_vector[i] <- if (normalize) {
      KT_calc_gini_norm(sampled_df$actual, sampled_df$weight, sampled_df$predicted)
    } else {
      KT_calc_gini(sampled_df$actual, sampled_df$weight, sampled_df$predicted)
    }
  }
  return(gini_vector)
}

# Plot Gini comparison
KT_plot_compare_gini <- function(n, actual, weight, base, challenger, normalize = FALSE) {
  base_gini <- KT_resample_gini(n, actual, weight, base, normalize)
  challenger_gini <- KT_resample_gini(n, actual, weight, challenger, normalize)
  challenger_win_rate <- mean(challenger_gini > base_gini)
  gini_df <- data.frame(Model = rep(c("Base", "Challenger"), each = n), Gini = c(base_gini, challenger_gini))
  ggplot(gini_df, aes(x = Gini, fill = Model)) +
    geom_density(alpha = 0.3) +
    ggtitle(glue("Gini Comparison | Challenger Win Rate: {scales::percent(challenger_win_rate)}"))
}

# Weighted R-squared
KT_weighted_Rsq <- function(actual, pred, weight) {
  residual_ss <- sum(((actual - pred)^2) * weight)
  total_ss <- sum(((actual - mean(actual))^2) * weight)
  return(1 - residual_ss / total_ss)
}

# Lift calculations
KT_calc_lift <- function(pred, actual, weight, nbin) {
  pred <- pred * (sum(actual) / sum(pred * weight))  # Rebase predictions
  lift_df <- data.frame(pred, actual, weight) %>%
    filter(weight > 0) %>%
    arrange(pred) %>%
    mutate(pred = pred * weight, bin = KT_create_equal_bin(weight, nbin)) %>%
    group_by(bin) %>%
    summarise(across(everything(), sum)) %>%
    mutate(actual = actual / weight, pred = pred / weight, AvE = actual / pred)
  return(lift_df)
}

# Double-lift calculations
KT_calc_dl <- function(actual, weight, base, challenger, nbin) {
  df <- data.frame(actual, weight, base, challenger) %>%
    filter(weight > 0) %>%
    mutate(model_ratio = base / challenger) %>%
    arrange(model_ratio) %>%
    mutate(bin = KT_create_equal_bin(weight, nbin)) %>%
    group_by(bin) %>%
    summarise(across(everything(), sum)) %>%
    mutate(actual = actual / weight, base = base / weight, challenger = challenger / weight, AvE = actual / base)
  return(df)
}
