# -------------------------
# Load libraries
# -------------------------
# Check if 'rstan' is installed, install if missing, then load
if (!requireNamespace("rstan", quietly = TRUE)) {
  install.packages("rstan", dependencies = TRUE)
}
# Load libraries
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# -------------------------
# 1. Load data
# -------------------------
# Example: replace with your actual data vectors
# Hard-coded dose vector
dose <- c(
  0.000000, 60.768000, 120.615273, 6.445091, 67.213091, 127.060364, 187.828364,
  0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
  0.000000, 13.810909, 73.658182, 134.426182, 194.273455, 20.256000, 80.103273,
  140.871273, 201.639273, 26.701091, 87.469091, 147.316364, 208.084364, 33.146182,
  93.914182, 154.682182, 214.529455, 100.359273, 161.127273, 220.974545, 40.512000,
  107.725091, 167.572364, 228.340364, 46.957091, 114.170182, 174.017455, 234.785455,
  53.402182
)

# Hard-coded response vector
response <- c(
  -0.007896428, 0.308991995, 0.795676998, 0.009178554, 0.477518552, 0.833293128,
  0.836721730, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000,
  0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.037092345, 0.486274322,
  1.000000000, 0.914326675, 0.067678572, 0.269712685, 0.733468553, 1.000000000,
  0.087855619, 0.309578124, 0.679342900, 0.857521709, 0.089436680, 0.297529801,
  0.753766669, 0.857291917, 0.527127213, 0.824219442, 0.922556208, 0.120712277,
  0.656446556, 0.700586086, 0.776056829, 0.143432026, 0.550213197, 0.774477511,
  1.000000000, 0.154565341
)

data_list <- list(
  N = length(dose),
  dose = dose,
  y = response
)

# -------------------------
# 2. Stan model code (heteroscedastic)
# -------------------------
stan_code <- "
data {
  int<lower=1> N;
  vector[N] dose;
  vector[N] y;
}
parameters {
  real<lower=0> alpha;
  real<lower=0> beta;
  real<lower=0,upper=1> gamma;
  real<lower=0> sigma0;        // base sigma
  real<lower=0> sigma_pow;     // exponent for dose dependence
}
model {
  vector[N] mu;
  vector[N] sigma;

  // Priors
  alpha ~ normal(100, 100);
  beta ~ normal(1, 5);
  gamma ~ beta(2, 2);
  sigma0 ~ normal(0, 0.1);
  sigma_pow ~ normal(0, 1);

  // Mean and dose-dependent variance
  for (i in 1:N) {
    mu[i] = gamma * (1 - exp(-pow(dose[i]/alpha, beta)));
    sigma[i] = sigma0 * pow(1 + dose[i], sigma_pow); // sigma increases with dose
  }

  // Likelihood
  y ~ normal(mu, sigma);
}
generated quantities {
  vector[N] y_pred;
  vector[N] sigma_pred;
  for (i in 1:N) {
    sigma_pred[i] = sigma0 * pow(1 + dose[i], sigma_pow);
    y_pred[i] = normal_rng(gamma * (1 - exp(-pow(dose[i]/alpha, beta))), sigma_pred[i]);
  }
}
"

# -------------------------
# 3. Fit the model
# -------------------------
fit <- stan(
  model_code = stan_code,
  data = data_list,
  iter = 4000,
  chains = 4,
  seed = 123
)

# -------------------------
# 4. Extract posterior predictive draws
# -------------------------
posterior <- rstan::extract(fit)
y_pred <- posterior$y_pred

quantiles <- apply(y_pred, 2, function(x) quantile(x, probs = c(0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95)))

# -------------------------
# 5. Plot with ggplot2
# -------------------------
plot_data <- data.frame(
  dose = dose,
  response = response,
  lower = quantiles[1,],
  LL = quantiles[2,],
  median = quantiles[3,],
  UU = quantiles[4,],
  upper = quantiles[5,]
)

# -------------------------
# 6. Plot in base R
# -------------------------
dose_grid <- seq(min(dose), max(dose), length.out = 200)

# Interpolate quantiles for smooth curve
lower_smooth <- approx(dose, quantiles[1,], xout = dose_grid)$y
LL_smooth <- approx(dose, quantiles[2,], xout = dose_grid)$y
ll_smooth <- approx(dose, quantiles[3,], xout = dose_grid)$y
median_smooth <- approx(dose, quantiles[4,], xout = dose_grid)$y
ul_smooth <- approx(dose, quantiles[5,], xout = dose_grid)$y
UL_smooth <- approx(dose, quantiles[6,], xout = dose_grid)$y
upper_smooth <- approx(dose, quantiles[7,], xout = dose_grid)$y

quartz(width = 205.5213/50, height = 139.6644/50)
# Base R plot
par(yaxs='r', xaxs='i', mar=c(4,4,1,1))
plot(dose, response, pch = 16, col = "black",
     xlab = "UV power (µW)",
     ylab = "Fluorescent labeling intensity (A.U.)",
     ylim = c(0, 1), las=1, axes=F)

axis(1)
axis(2, at=c(0, 0.5, 1), las=1)

# Add smooth ribbon
polygon(c(dose_grid, rev(dose_grid)),
        c(upper_smooth, rev(lower_smooth)),
        col = rgb(0.68, 0.85, 0.90, alpha = 0.5), border = 'black', lty=2)

polygon(c(dose_grid, rev(dose_grid)),
        c(UL_smooth, rev(LL_smooth)),
        col = rgb(0.68, 0.85, 0.90, alpha = 0.5), border = 'black', lty=3)

polygon(c(dose_grid, rev(dose_grid)),
        c(ul_smooth, rev(ll_smooth)),
        col = rgb(0.68, 0.85, 0.90, alpha = 0.5), border = 'black', lty=3)

# Add vertical error bar segments
segments(
  x0 = dose,
  y0 = response - merged$Std_norm,
  x1 = dose,
  y1 = response + merged$Std_norm,
  col = "black",
  lwd = 1.5
)


points(dose, response, pch = 21, bg='white', cex=1.5)

# Add median line
lines(dose_grid, median_smooth, col = "red", lwd = 2)
box()

quartz.save(file="./figures/pdf/intensity_SuppFigure07_panelE.pdf", type='pdf')
