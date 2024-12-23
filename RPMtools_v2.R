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

# Options
options(scipen = 999)

# Utility: Create Equal Bins
KT_create_equal_bin <- function(weight, nbin) {
  if (any(is.na(weight)) || any(weight < 0)) stop("Invalid weights: Weights must be non-negative and not NA.")
  if (length(weight) == 0) stop("Weight vector is empty.")
  if (nbin <= 0) stop("Number of bins must be a positive integer.")
  
  cumulative_sum <- cumsum(weight)
  bins <- cut(cumulative_sum, breaks = nbin, labels = FALSE)
  if (any(is.na(bins))) warning("Some bins are empty due to insufficient data.")
  return(bins)
}

# Gini Calculation
KT_calc_gini <- function(actual, weight, predicted) {
  if (length(actual) != length(weight) || length(actual) != length(predicted)) {
    stop("Input vectors 'actual', 'weight', and 'predicted' must have the same length.")
  }
  if (any(is.na(actual)) || any(is.na(weight)) || any(is.na(predicted))) {
    stop("Inputs contain NA values. Please remove or impute them.")
  }
  if (all(weight == 0)) stop("All weights are zero. Cannot compute Gini.")
  
  df <- data.frame(actual = as.numeric(actual), weight = as.numeric(weight), predicted)
  sorted_idx <- order(df$predicted)
  w_s <- df$weight[sorted_idx]
  a_s <- df$actual[sorted_idx]
  a_c <- cumsum(a_s * w_s)
  w_c <- cumsum(w_s)
  gini <- 1 - 2 * pracma::trapz(w_c / max(w_c), a_c / max(a_c))
  return(gini)
}

# Normalized Gini
KT_calc_gini_norm <- function(actual, weight, predicted) {
  return(KT_calc_gini(actual, weight, predicted) / KT_calc_gini(actual, weight, actual))
}

# Resample Gini
KT_resample_gini <- function(n, actual, weight, predicted, normalize = FALSE) {
  if (n <= 0) stop("Number of resamples must be a positive integer.")
  if (any(weight < 0)) warning("Negative weights detected; results may be incorrect.")
  
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

# Gini Plotting
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
  if (length(actual) != length(pred) || length(actual) != length(weight)) {
    stop("Input vectors 'actual', 'pred', and 'weight' must have the same length.")
  }
  if (any(is.na(actual)) || any(is.na(pred)) || any(is.na(weight))) {
    stop("Inputs contain NA values. Please remove or impute them.")
  }
  if (all(weight == 0)) stop("All weights are zero. Cannot compute weighted R-squared.")
  
  residual_ss <- sum(((actual - pred)^2) * weight)
  total_ss <- sum(((actual - mean(actual))^2) * weight)
  r_squared <- 1 - residual_ss / total_ss
  return(r_squared)
}

# Lift Calculation
KT_calc_lift <- function(pred, actual, weight, nbin) {
  if (length(pred) != length(actual) || length(actual) != length(weight)) {
    stop("Input vectors 'pred', 'actual', and 'weight' must have the same length.")
  }
  if (nbin <= 0) stop("Number of bins must be a positive integer.")
  
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

# Double Lift Calculation
KT_calc_dl <- function(actual, weight, base, challenger, nbin) {
  if (length(actual) != length(weight) || length(actual) != length(base) || length(actual) != length(challenger)) {
    stop("Input vectors 'actual', 'weight', 'base', and 'challenger' must have the same length.")
  }
  if (nbin <= 0) stop("Number of bins must be a positive integer.")
  
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
################################ AvE ###################################

# Function: Calculate Actual vs Expected (AvE)
KT_calc_ave <- function(ft, actual, pred, challenger, weight, rebase = TRUE) {
  if (missing(ft)) stop("Feature ('ft') is missing.")
  if (missing(actual) || missing(pred) || missing(weight)) stop("One or more required inputs ('actual', 'pred', 'weight') are missing.")
  if (any(is.na(c(ft, actual, pred, weight)))) stop("Inputs contain NA values. Please remove or impute them.")
  
  if (missing(challenger)) challenger <- pred

  if (rebase) {
    pred <- pred * (sum(actual) / sum(pred * weight))
    challenger <- challenger * (sum(actual) / sum(challenger * weight))
  }

  df <- data.frame(ft, actual, pred, challenger, weight)
  df <- df %>%
    mutate(across(c(pred, challenger), ~ .x * weight))  # Apply weights
  
  overall <- df %>%
    summarise(across(c(actual, pred, challenger, weight), sum)) %>%
    mutate(
      actual_overall_avg = actual / weight,
      pred_overall_avg = pred / weight,
      challenger_overall_avg = challenger / weight
    )
  
  result <- df %>%
    group_by(ft) %>%
    summarise(across(c(actual, pred, challenger, weight), sum)) %>%
    mutate(
      actual = actual / weight,
      pred = pred / weight,
      challenger = challenger / weight,
      ave = actual / pred,
      challenger_ave = actual / challenger,
      actual_overall_avg = overall$actual_overall_avg,
      pred_overall_avg = overall$pred_overall_avg
    )
  
  return(result)
}

# Function: Random Fold AvE Consistency
KT_calc_ave_consistency_random_fold <- function(ft, actual, pred, weight, challenger, nfold = 5, plot_scale = 5000) {
  if (missing(ft)) stop("Feature ('ft') is missing.")
  if (nfold <= 0) stop("Number of folds ('nfold') must be a positive integer.")
  
  folds <- KT_create_fold_idx(data.frame(ft), nfold)
  folds[["Full"]] <- unlist(folds) %>% as.vector()
  AvE_df_list <- list()

  for (fold in names(folds)) {
    fold_data <- folds[[fold]]
    AvE_df_list[[fold]] <- KT_calc_ave(ft = ft[fold_data], actual = actual[fold_data], pred = pred[fold_data], weight = weight[fold_data]) %>%
      mutate(sample = fold)
  }

  ave_df <- rbindlist(AvE_df_list) %>%
    mutate(
      sample = factor(sample, levels = KT_dym_sort(unique(sample))),
      bar_group = ifelse(grepl("fold", sample), "fold", "full")
    )
  
  p <- ggplotly(
    ggplot(ave_df, aes(x = ft, group = sample, fill = bar_group)) +
      geom_hline(yintercept = 1, color = "#39ff14") +
      geom_line(aes(y = ave, color = sample)) +
      geom_bar(aes(y = weight / plot_scale), stat = "identity", alpha = 0.4, position = "dodge") +
      scale_fill_manual(values = c("fold" = "grey", "full" = "orange")) +
      scale_y_continuous(name = "Actual/Expected", sec.axis = sec_axis(~ . * plot_scale, name = "weight")) +
      theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 0.9)) +
      ggtitle("AvE Consistency Across Random Folds")
  )
  
  return(list(ave_df = ave_df, ave_plot = p))
}

# Function: Resample AvE
KT_resample_ave <- function(n, ft, actual, pred, challenger, weight) {
  if (n <= 0) stop("Number of resamples ('n') must be a positive integer.")
  if (missing(challenger)) challenger <- pred

  ave_sim <- list()
  df <- data.frame(ft, actual, pred, challenger, weight)
  
  main_ave <- KT_calc_ave(ft = df$ft, actual = df$actual, pred = df$pred, challenger = df$challenger, weight = df$weight)
  main_ave$sample <- "main"
  ave_sim[["iter_0"]] <- main_ave

  for (x in seq_len(n)) {
    set.seed(x)
    sampled_df <- df %>% sample_frac(size = 0.3, replace = FALSE)
    ave_sim[[glue("iter_{x}")]] <- KT_calc_ave(sampled_df$ft, sampled_df$actual, sampled_df$pred, sampled_df$challenger, sampled_df$weight) %>%
      mutate(sample = x)
  }

  variables <- list()
  for (var in c("actual", "pred", "ave", "challenger", "challenger_ave")) {
    variables[[var]] <- rbindlist(ave_sim) %>%
      select(ft, !!as.name(var), sample) %>%
      pivot_wider(names_from = sample, values_from = !!as.name(var)) %>%
      rowwise() %>%
      mutate(lb = quantile(c_across(2:(n + 1)), 0.05, na.rm = TRUE),
             ub = quantile(c_across(2:(n + 1)), 0.95, na.rm = TRUE)) %>%
      select(ft, main, lb, ub) %>%
      mutate(variable = var)
  }

  ave_df <- data.frame(rbindlist(variables), weight = main_ave$weight)
  return(list(ave_df = ave_df, main_ave = main_ave))
}
