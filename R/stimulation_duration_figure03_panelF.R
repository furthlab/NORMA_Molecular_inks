# ==============================================================================
# Figure 3f | Stimulation duration and labeling intensity
# ==============================================================================

library(sciplot)


# ---- Parameters --------------------------------------------------------------

data_dir <- "./data/csv/stimulation_duration"

output_file <- "./figures/pdf/figure03_panelF_stimulation_duration.pdf"

duration_levels <- c(
  "background", "1ms", "5ms", "15ms", "31ms",
  "63ms", "125ms", "250ms", "500ms", "1000ms"
)

excluded_durations_ms <- c(1, 5)
confidence_level <- 0.999

background_color <- "gray"
stimulation_color <- "darkred"

figure_height <- 3.835294
figure_width  <- (433.9209 / 198.3304) * figure_height


# ---- Load experimental data --------------------------------------------------

files <- list.files(
  data_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)

stopifnot(length(files) > 0)

dat <- do.call(
  rbind,
  lapply(files, function(file) {
    df <- read.csv(file)
    df$dur <- sub("\\..*$", "", basename(file))
    df
  })
)

# Define stimulation duration order
dat$dur <- factor(dat$dur, levels = duration_levels)

# Remove missing or non-finite measurements
dat <- subset(dat, !is.na(dur) & is.finite(Mean))

# Calculate mean background intensity
background_mean <- mean(dat$Mean[dat$dur == "background"])


# ---- Calculate descriptive statistics ----------------------------------------

group_stats <- do.call(
  rbind,
  lapply(duration_levels, function(duration) {

    values <- dat$Mean[dat$dur == duration]

    n <- length(values)

    data.frame(
      Duration = duration,
      N        = n,
      Mean     = mean(values),
      Median   = median(values),
      SD       = sd(values),
      SEM      = sd(values) / sqrt(n)
    )
  })
)

group_stats$Included_in_regression <-
  !group_stats$Duration %in% c(
    "background",
    paste0(excluded_durations_ms, "ms")
  )


# ---- Prepare regression data -------------------------------------------------

# Exclude background measurements
dat_stim <- subset(dat, dur != "background")

# Convert stimulation duration to milliseconds
dat_stim$dur_ms <- as.numeric(
  sub("ms", "", as.character(dat_stim$dur))
)

# Select measurements contributing to regression
dat_regression <- subset(
  dat_stim,
  !dur_ms %in% excluded_durations_ms
)

# Calculate arithmetic mean intensity for each duration
means <- aggregate(
  Mean ~ dur_ms,
  data = dat_regression,
  FUN = mean
)

# Verify that all group means are positive before log transformation
stopifnot(
  all(means$Mean > 0),
  all(means$dur_ms > 0)
)

# Fit linear regression to log10-transformed group means
fit <- lm(
  log10(Mean) ~ log10(dur_ms),
  data = means
)

fit_summary <- summary(fit)

coefs <- coef(fit)
r2 <- fit_summary$r.squared

regression_label <- sprintf(
  "log10(Labeling) = %.1f + %.1f log10(Duration)\nR² = %.2f",
  coefs[1],
  coefs[2],
  r2
)


# ---- Calculate regression confidence interval --------------------------------

prediction_x <- seq(0, 4, length.out = 100)

confidence_interval <- predict(
  fit,
  newdata = data.frame(dur_ms = 10^prediction_x),
  interval = "confidence",
  level = confidence_level
)


# ---- Report descriptive statistics -------------------------------------------

cat("\n")
cat("============================================================\n")
cat("FIGURE 3F: DESCRIPTIVE STATISTICS\n")
cat("============================================================\n\n")

print(
  group_stats,
  row.names = FALSE,
  digits = 5
)

cat("\nTotal number of measurements:", nrow(dat), "\n")

cat(
  "Background measurements:",
  sum(dat$dur == "background"), "\n"
)

cat(
  "Stimulation measurements:",
  nrow(dat_stim), "\n"
)

cat(
  "Number of stimulation durations:",
  length(unique(dat_stim$dur_ms)), "\n"
)

cat(
  "Mean background intensity:",
  background_mean, "\n"
)


# ---- Report regression sample sizes ------------------------------------------

cat("\n")
cat("============================================================\n")
cat("REGRESSION SAMPLE SIZES\n")
cat("============================================================\n\n")

regression_stats <- subset(
  group_stats,
  Included_in_regression
)

print(
  regression_stats,
  row.names = FALSE,
  digits = 5
)

cat(
  "\nIncluded stimulation durations (ms):",
  paste(means$dur_ms, collapse = ", "), "\n"
)

cat(
  "Excluded stimulation durations (ms):",
  paste(excluded_durations_ms, collapse = ", "), "\n"
)

cat(
  "Number of underlying measurements contributing to group means:",
  nrow(dat_regression), "\n"
)

cat(
  "Number of excluded stimulation measurements:",
  nrow(dat_stim) - nrow(dat_regression), "\n"
)

cat(
  "Number of group means fitted in regression:",
  nobs(fit), "\n"
)

cat(
  "Number of fitted regression parameters:",
  length(coef(fit)), "\n"
)

cat(
  "Residual degrees of freedom:",
  df.residual(fit), "\n"
)


# ---- Report regression statistics --------------------------------------------

cat("\n")
cat("============================================================\n")
cat("LOG-LOG LINEAR REGRESSION STATISTICS\n")
cat("============================================================\n\n")

print(fit_summary)

# Regression coefficients and their standard errors
coefficient_stats <- as.data.frame(fit_summary$coefficients)

colnames(coefficient_stats) <- c(
  "Estimate",
  "SE",
  "t_value",
  "p_value"
)

# Calculate parameter confidence intervals
parameter_ci <- confint(
  fit,
  level = confidence_level
)

coefficient_stats$CI_lower <- parameter_ci[, 1]
coefficient_stats$CI_upper <- parameter_ci[, 2]

cat("\nRegression coefficients:\n\n")

print(coefficient_stats, digits = 6)

cat(
  "\nR-squared:",
  fit_summary$r.squared, "\n"
)

cat(
  "Adjusted R-squared:",
  fit_summary$adj.r.squared, "\n"
)

cat(
  "Residual standard error:",
  fit_summary$sigma, "\n"
)

cat(
  "Confidence interval level (%):",
  confidence_level * 100, "\n"
)

# Overall regression significance
f_stat <- fit_summary$fstatistic

model_p <- pf(
  f_stat[1],
  f_stat[2],
  f_stat[3],
  lower.tail = FALSE
)

cat("F-statistic:", unname(f_stat[1]), "\n")
cat("Numerator degrees of freedom:", unname(f_stat[2]), "\n")
cat("Denominator degrees of freedom:", unname(f_stat[3]), "\n")
cat("Overall regression p-value:", model_p, "\n")


# ---- Report individual regression observations -------------------------------

cat("\n")
cat("============================================================\n")
cat("INDIVIDUAL REGRESSION OBSERVATIONS\n")
cat("============================================================\n\n")

regression_observations <- data.frame(
  Duration_ms = means$dur_ms,
  Mean_intensity = means$Mean,
  Log10_duration = log10(means$dur_ms),
  Log10_intensity = log10(means$Mean),
  Fitted_log10_intensity = fitted(fit),
  Residual = residuals(fit)
)

print(
  regression_observations,
  row.names = FALSE,
  digits = 6
)


# ---- Prepare logarithmic axes ------------------------------------------------

log_minor_ticks <- function(exponents) {
  log10(as.vector(outer(1:9, 10^exponents)))
}

x_minor_ticks <- log_minor_ticks(0:3)
y_minor_ticks <- log_minor_ticks(1:3)

x_major_ticks <- 0:4
y_major_ticks <- 1:4

x_labels <- parse(text = paste0("10^", x_major_ticks))
y_labels <- parse(text = paste0("10^", y_major_ticks))


# ---- Initialize figure -------------------------------------------------------

quartz(
  width = figure_width,
  height = figure_height
)

par(
  mfrow = c(1, 2),
  mar = c(4, 4, 1, 1),
  yaxs = "i"
)


# ---- Left panel: Labeling intensity by stimulation duration -------------------

bargraph.CI(
  dat$dur,
  dat$Mean,
  ci.fun = function(x) {
    c(
      mean(x) - sd(x) / sqrt(length(x)),
      mean(x) + sd(x) / sqrt(length(x))
    )
  },
  ylab = "Labeling intensity",
  xlab = "Stimulation duration (ms)",
  col = c(background_color, rep(stimulation_color, 9)),
  las = 1,
  ylim = c(0, 1200)
)

# Background reference
abline(
  h = background_mean,
  lty = 2,
  lwd = 3
)


# ---- Right panel: Log-log regression -----------------------------------------

plot(
  log10(means$dur_ms),
  log10(means$Mean),
  type = "n",
  axes = FALSE,
  xlim = c(0, 4),
  ylim = c(1, 4),
  xlab = "Stimulation duration (ms)",
  ylab = "Labeling intensity"
)

# Minor logarithmic ticks
axis(
  side = 1,
  at = x_minor_ticks,
  labels = FALSE,
  tck = -0.018
)

axis(
  side = 2,
  at = y_minor_ticks,
  labels = FALSE,
  tck = -0.018
)

# Major logarithmic ticks
axis(
  side = 1,
  at = x_major_ticks,
  labels = x_labels
)

axis(
  side = 2,
  at = y_major_ticks,
  labels = y_labels,
  las = 1
)

# Fitted regression line
abline(fit, lwd = 2)

# 99.9% confidence interval for the fitted mean
matlines(
  prediction_x,
  confidence_interval[, c("lwr", "upr")],
  col = "blue",
  lty = 2
)

# Mean labeling intensity at each included duration
points(
  log10(means$dur_ms),
  log10(means$Mean),
  pch = 21,
  cex = 2,
  bg = "white"
)

# Regression equation and coefficient of determination
text(
  x = 0,
  y = 4,
  labels = regression_label,
  adj = c(0, 2),
  cex = 0.8
)


# ---- Export figure -----------------------------------------------------------

dir.create(
  dirname(output_file),
  recursive = TRUE,
  showWarnings = FALSE
)

quartz.save(
  file = output_file,
  type = "pdf"
)

cat("\nFigure saved to:", output_file, "\n")
