# ==============================================================================
# Figure 3d | Spatial resolution of Gaussian stripe photopatterning
# ==============================================================================

# ---- Parameters --------------------------------------------------------------
# https://colab.research.google.com/drive/1w-cs25UL9ELCUWsfXOLTDFN6q7tBKCm1#scrollTo=lF2JHMbyQioh

# Theoretical Gaussian stripe pattern
pattern_pixel_size <- 0.1                  # µm/pixel
pattern_width_um   <- 400                  # µm
stripe_sigmas_um   <- c(0.2, 0.4, 0.8, 1.6, 3.2, 6.4)
spacing_factor     <- 6

# Experimental imaging
pixel_size <- 0.1624                       # µm/pixel

# Spatial alignment parameters
theoretical_scale  <- 2.05
theoretical_offset <- -140
labeling_scale     <- 0.99
labeling_offset    <- -240

# Intensity scaling and offsets
theoretical_gain     <- 1000
theoretical_baseline <- 1400
labeling_baseline    <- 1350


# ---- Generate theoretical Gaussian stripe pattern ---------------------------

# Calculate proportional spacing between adjacent stripes
stripe_spacing_um <- spacing_factor *
  (head(stripe_sigmas_um, -1) + tail(stripe_sigmas_um, -1)) / 2

centers_um <- c(0, cumsum(stripe_spacing_um))

# Center the stripe pattern horizontally
centers_um <- centers_um - mean(range(centers_um)) +
  pattern_width_um / 2

# Generate spatial coordinates
x_um <- seq(
  from = 0,
  by = pattern_pixel_size,
  length.out = as.integer(pattern_width_um / pattern_pixel_size)
)

# Sum Gaussian intensity profiles
intensity <- Reduce(
  `+`,
  Map(
    function(center, sigma) {
      exp(-0.5 * ((x_um - center) / sigma)^2)
    },
    centers_um,
    stripe_sigmas_um
  )
)

# Normalize theoretical intensity to [0, 1]
intensity <- intensity / max(intensity)


# ---- Prepare theoretical intensity profile ----------------------------------

# Extract region containing the Gaussian stripe array
theoretical_idx <- 1585:2550
theoretical_y <- intensity[theoretical_idx]

theoretical_x <- seq(
  from = 1,
  to = length(theoretical_y) * theoretical_scale,
  length.out = length(theoretical_y)
) + theoretical_offset

# Identify theoretical stripe centers within the extracted profile
peak_idx <- vapply(
  centers_um,
  function(center) {
    which.min(abs(x_um[theoretical_idx] - center))
  },
  integer(1)
)

peak_positions <- theoretical_x[peak_idx]


# ---- Load experimental intensity profiles -----------------------------------

# UV illumination profile
uv_data <- read.csv(
  "./data/csv/UV_stim_gaussian_stripe_array.csv"
)

uv_median <- apply(uv_data[, -1], 1, median)

# Fluorescent labeling profile
labeling_data <- read.csv(
  "./data/csv/labeled_gaussian_stripe_array.csv"
)

# Median and standard deviation across 39 technical replicates
labeling_median <- apply(labeling_data[, -1], 1, median)
labeling_sd     <- apply(labeling_data[, -1], 1, sd)

# Spatially align the labeling profile
labeling_x <- seq(
  from = 1,
  to = length(labeling_median) * labeling_scale,
  length.out = length(labeling_median)
) + labeling_offset

# Apply intensity offset
labeling_y <- labeling_median - labeling_baseline

# Calculate median ± SD
labeling_upper <- labeling_y + labeling_sd
labeling_lower <- labeling_y - labeling_sd


# ---- Prepare spatial axis ----------------------------------------------------

micron_ticks <- seq(
  from = 0,
  to = floor(1967 * pixel_size / 50) * 50,
  by = 50
)

pixel_ticks <- micron_ticks / pixel_size


# ---- Plot -------------------------------------------------------------------

quartz(
  width = 208.2416 / 30,
  height = 85.084 / 30
)

par(mar = c(4, 3, 0, 0))

# Experimental UV illumination profile
plot(
  head(uv_median, -1),
  type = "l",
  ylim = c(-1100, 2500),
  xlim = c(50, 1650),
  col = "purple",
  lwd = 2,
  axes = FALSE,
  ylab = "",
  xlab = "Distance (µm)"
)

# Theoretical Gaussian stripe profile
lines(
  theoretical_x,
  theoretical_y * theoretical_gain + theoretical_baseline,
  lwd = 2
)

# Fluorescent labeling variability (median ± SD)
polygon(
  c(labeling_x, rev(labeling_x)),
  c(labeling_upper, rev(labeling_lower)),
  col = "#FA825F80",
  border = "#D0416F"
)

# Median fluorescent labeling profile
lines(
  labeling_x,
  labeling_y,
  col = "black",
  lwd = 4
)

lines(
  labeling_x,
  labeling_y,
  col = "orange",
  lwd = 2
)

# Gaussian stripe annotations
axis(
  side = 1,
  at = peak_positions,
  labels = c(0, 6, 12, 24, 48, 96)
)

# Physical distance axis
axis(
  side = 1,
  at = pixel_ticks,
  labels = micron_ticks,
  line = 2
)


# ---- Export figure -----------------------------------------------------------

dir.create(
  "./figures/pdf",
  recursive = TRUE,
  showWarnings = FALSE
)

quartz.save(
  file = "./figures/pdf/figure03_panelD_spatial_resolution.pdf",
  type = "pdf"
)
