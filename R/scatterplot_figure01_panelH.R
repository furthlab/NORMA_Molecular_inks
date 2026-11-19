library(MASS)
library(stringi)

# ============================================================
# SETTINGS
# ============================================================

data_dir <- "./data/microscope/multiplex"

area_cutoff <- 1000

xlim_plot <- c(6.5, 12)
ylim_plot <- c(6.5, 12)

kde_resolution <- 300

point_size <- 1.5


# ============================================================
# FIND CSV FILES
# ============================================================

files <- list.files(
  data_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)

print(basename(files))


# ============================================================
# READ EACH FILE
# ============================================================

read_file <- function(f) {

  d <- read.csv(f)

  # Extract channel number from Label
  channel <- stringi::stri_split_fixed(
    d$Label,
    ":",
    simplify = TRUE
  )[, 3]


  # ----------------------------------------------------------
  # Determine experimental condition from filename
  # ----------------------------------------------------------

  if (grepl("MAX_U2OS_MCDIBOD_DAMA", basename(f))) {

    condition <- "Both"

  } else if (grepl("MAX_20060915_U2OS_MCDIBOD_488", basename(f))) {

    condition <- "NORMA315"

  } else if (grepl("MAX_20060915_U2OS_NORMA345", basename(f))) {

    condition <- "NORMA345"

  } else {

    message("Skipping unknown condition: ", basename(f))
    return(NULL)
  }


  # ----------------------------------------------------------
  # Find rows belonging to channel 1 and channel 2
  # ----------------------------------------------------------

  i1 <- which(channel == "1")
  i2 <- which(channel == "2")


  if (length(i1) == 0 || length(i2) == 0) {

    message(
      "Skipping file because a channel is missing: ",
      basename(f)
    )

    return(NULL)
  }


  # ----------------------------------------------------------
  # Pair channel 1 and channel 2 objects
  # ----------------------------------------------------------

  n <- min(
    length(i1),
    length(i2)
  )

  i1 <- i1[seq_len(n)]
  i2 <- i2[seq_len(n)]


  # ----------------------------------------------------------
  # Area filter
  # ----------------------------------------------------------

  keep <- d$Area[i1] > area_cutoff

  keep[is.na(keep)] <- FALSE


  if (!any(keep)) {

    message(
      "Skipping file because no objects pass area cutoff: ",
      basename(f)
    )

    return(NULL)
  }


  # ----------------------------------------------------------
  # Extract intensities
  # ----------------------------------------------------------

  x <- log2(
    d$Mean[i1][keep]
  )

  y <- log2(
    d$Mean[i2][keep]
  )

  area <- d$Area[i1][keep]


  # Remove NA / Inf / -Inf
  valid <- (
    is.finite(x) &
      is.finite(y) &
      is.finite(area)
  )

  x <- x[valid]
  y <- y[valid]
  area <- area[valid]


  if (length(x) == 0) {

    message(
      "Skipping file because no finite values remain: ",
      basename(f)
    )

    return(NULL)
  }


  # ----------------------------------------------------------
  # Return paired data
  # ----------------------------------------------------------

  data.frame(
    x = x,
    y = y,
    Area = area,
    file = rep(basename(f), length(x)),
    condition = rep(condition, length(x))
  )
}


# ============================================================
# COMBINE ALL FILES
# ============================================================

data_list <- lapply(
  files,
  read_file
)

data_list <- data_list[
  !sapply(data_list, is.null)
]

data <- do.call(
  rbind,
  data_list
)


# ============================================================
# CHECK NUMBER OF POINTS
# ============================================================

cat("\nPoints per condition:\n")

print(
  table(data$condition)
)

cat("\nPoints per file:\n")

print(
  table(data$file)
)


# ============================================================
# SAVE RAW VALUES
# ============================================================

data$x_raw <- data$x
data$y_raw <- data$y


# ============================================================
# SHOW RAW REPLICATE MEDIANS
# ============================================================

cat("\nRaw replicate medians, channel X:\n")

print(
  aggregate(
    x_raw ~ condition + file,
    data = data,
    FUN = median
  )
)

cat("\nRaw replicate medians, channel Y:\n")

print(
  aggregate(
    y_raw ~ condition + file,
    data = data,
    FUN = median
  )
)


# ============================================================
# NORMALIZE REPLICATE EFFECTS
#
# Normalization happens ONLY within each experimental condition.
#
# For each condition:
#
#   1. Calculate median X and Y for each replicate.
#   2. Calculate the median of the replicate medians.
#   3. Shift every replicate so its median matches that target.
#
# Because intensities are already log2 transformed, this
# corresponds to multiplicative normalization in raw intensity.
#
# Importantly:
#
# NORMA315, NORMA345 and Both are NOT normalized to each other.
# Their differences in overall position are retained.
# ============================================================

data$x_norm <- NA_real_
data$y_norm <- NA_real_


for (cond in unique(data$condition)) {

  # Data belonging to this condition
  ii_cond <- data$condition == cond

  condition_files <- unique(
    data$file[ii_cond]
  )


  # ----------------------------------------------------------
  # Calculate replicate medians
  # ----------------------------------------------------------

  replicate_x_medians <- sapply(
    condition_files,
    function(f) {

      median(
        data$x_raw[
          data$condition == cond &
            data$file == f
        ],
        na.rm = TRUE
      )

    }
  )


  replicate_y_medians <- sapply(
    condition_files,
    function(f) {

      median(
        data$y_raw[
          data$condition == cond &
            data$file == f
        ],
        na.rm = TRUE
      )

    }
  )


  # ----------------------------------------------------------
  # Target = median of replicate medians
  #
  # This gives each biological/image replicate equal weight,
  # regardless of how many cells/objects are in that file.
  # ----------------------------------------------------------

  target_x <- median(
    replicate_x_medians,
    na.rm = TRUE
  )

  target_y <- median(
    replicate_y_medians,
    na.rm = TRUE
  )


  cat(
    "\nCondition:",
    cond,
    "\n"
  )

  cat(
    "Target X median:",
    target_x,
    "\n"
  )

  cat(
    "Target Y median:",
    target_y,
    "\n"
  )


  # ----------------------------------------------------------
  # Shift each replicate to condition target
  # ----------------------------------------------------------

  for (f in condition_files) {

    ii <- (
      data$condition == cond &
        data$file == f
    )


    replicate_x <- median(
      data$x_raw[ii],
      na.rm = TRUE
    )

    replicate_y <- median(
      data$y_raw[ii],
      na.rm = TRUE
    )


    shift_x <- target_x - replicate_x
    shift_y <- target_y - replicate_y


    cat(
      basename(f),
      "\n  X shift:",
      round(shift_x, 3),
      "\n  Y shift:",
      round(shift_y, 3),
      "\n"
    )


    data$x_norm[ii] <-
      data$x_raw[ii] +
      shift_x

    data$y_norm[ii] <-
      data$y_raw[ii] +
      shift_y
  }
}


# ============================================================
# USE NORMALIZED VALUES FROM HERE ON
# ============================================================

data$x <- data$x_norm
data$y <- data$y_norm


# ============================================================
# CHECK NORMALIZATION
# ============================================================

cat("\nNormalized replicate medians, X:\n")

print(
  aggregate(
    x ~ condition + file,
    data = data,
    FUN = median
  )
)

cat("\nNormalized replicate medians, Y:\n")

print(
  aggregate(
    y ~ condition + file,
    data = data,
    FUN = median
  )
)


# ============================================================
# OPTIONAL DIAGNOSTIC:
# REPlicate distributions AFTER NORMALIZATION
# ============================================================

condition_colors <- c(
  "Both" = "grey40",
  "NORMA315" = "green",
  "NORMA345" = "magenta"
)

par(
  mfrow = c(1, 3),
  mar = c(4, 4, 2, 1)
)

for (cond in c(
  "Both",
  "NORMA315",
  "NORMA345"
)) {

  d <- data[
    data$condition == cond,
  ]

  ff <- unique(d$file)

  file_index <- match(
    d$file,
    ff
  )

  plot(
    d$x,
    d$y,
    pch = 16,
    cex = 0.35,
    col = file_index,
    asp = 1,
    las = 1,
    xlim = xlim_plot,
    ylim = ylim_plot,
    xlab = "NORMA-315",
    ylab = "NORMA-345",
    main = cond
  )

  legend(
    "topleft",
    legend = paste(
      "Replicate",
      seq_along(ff)
    ),
    col = seq_along(ff),
    pch = 16,
    bty = "n",
    cex = 0.8
  )
}

par(
  mfrow = c(1, 1)
)


# ============================================================
# CALCULATE KDE
#
# KDE is calculated independently for each condition using
# NORMALIZED x and y values.
# ============================================================

add_kde <- function(d) {

  # ----------------------------------------------------------
  # KDE limits based on actual data range
  # ----------------------------------------------------------

  xr <- range(
    d$x,
    finite = TRUE
  )

  yr <- range(
    d$y,
    finite = TRUE
  )


  # Add margin
  padx <- max(
    diff(xr) * 0.05,
    0.25
  )

  pady <- max(
    diff(yr) * 0.05,
    0.25
  )


  kde_limits <- c(
    xr[1] - padx,
    xr[2] + padx,
    yr[1] - pady,
    yr[2] + pady
  )


  # ----------------------------------------------------------
  # Calculate 2D KDE
  # ----------------------------------------------------------

  kd <- MASS::kde2d(
    d$x,
    d$y,
    n = kde_resolution,
    lims = kde_limits
  )


  # ----------------------------------------------------------
  # Find nearest KDE grid coordinate for each point
  # ----------------------------------------------------------

  ix <- findInterval(
    d$x,
    kd$x
  )

  iy <- findInterval(
    d$y,
    kd$y
  )


  ix <- pmax(
    1,
    pmin(
      ix,
      length(kd$x)
    )
  )

  iy <- pmax(
    1,
    pmin(
      iy,
      length(kd$y)
    )
  )


  # ----------------------------------------------------------
  # Extract density at each point
  # ----------------------------------------------------------

  d$density <- kd$z[
    cbind(ix, iy)
  ]


  return(d)
}


# ============================================================
# APPLY KDE SEPARATELY TO EACH CONDITION
# ============================================================

data <- do.call(
  rbind,
  lapply(
    split(
      data,
      data$condition
    ),
    add_kde
  )
)


# ============================================================
# CHECK KDE
# ============================================================

cat("\nNumber of NA KDE values:\n")

print(
  tapply(
    is.na(data$density),
    data$condition,
    sum
  )
)


cat("\nMean KDE density:\n")

print(
  tapply(
    data$density,
    data$condition,
    mean
  )
)


# ============================================================
# NORMALIZE KDE DENSITY WITHIN EACH CONDITION
# ============================================================

data$density_scaled <- NA_real_


for (cond in unique(data$condition)) {

  ii <- data$condition == cond

  z <- data$density[ii]

  zmin <- min(
    z,
    na.rm = TRUE
  )

  zmax <- max(
    z,
    na.rm = TRUE
  )


  if (zmax > zmin) {

    data$density_scaled[ii] <-
      (z - zmin) /
      (zmax - zmin)

  } else {

    data$density_scaled[ii] <- 1
  }
}


# ============================================================
# COLOR PALETTES
# ============================================================

palettes <- list(

  # Both NORMA-315 + NORMA-345 stimulated
  Both = colorRampPalette(
    c(
      "black",
      "white"
    )
  )(256),

  # NORMA-315 only
  NORMA315 = colorRampPalette(
    c(
      "darkgreen",
      "green"
    )
  )(256),

  # NORMA-345 only
  NORMA345 = colorRampPalette(
    c(
      "#8B008B",
      "magenta"
    )
  )(256)

)


# ============================================================
# FINAL PLOT
# ============================================================


quartz(height = 77.4092/20, width = 116.8296/20)

par(
  mar = c(6, 12, 1, 1)
)

plot(
  NA,
  xlim = xlim_plot,
  ylim = ylim_plot,
  asp = 1,
  las = 1,
  xlab = "NORMA-315",
  ylab = "NORMA-345",
  xaxt = "n",
  yaxt = "n"
)

# Major ticks: powers of 2
ticks <- 7:12

axis(
  1,
  at = ticks,
  labels = parse(text = paste0("2^", ticks))
)

axis(
  2,
  at = ticks,
  labels = parse(text = paste0("2^", ticks)),
  las = 1
)

# Minor ticks: linear subdivisions between each 2^n and 2^(n+1)
minor <- unlist(lapply(ticks[-length(ticks)], function(k) {
  log2(2^k * c(1.25, 1.5, 1.75))
}))

axis(
  1,
  at = minor,
  labels = FALSE,
  tck = -0.01
)

axis(
  2,
  at = minor,
  labels = FALSE,
  tck = -0.01
)

# ============================================================
# DRAW CONDITIONS
# ============================================================

for (cond in c(
  "Both",
  "NORMA315",
  "NORMA345"
)) {

  d <- data[
    data$condition == cond,
  ]


  if (nrow(d) == 0)
    next


  # Convert normalized KDE density to palette index
  col_idx <- round(
    d$density_scaled * 255
  ) + 1


  col_idx <- pmax(
    1,
    pmin(
      256,
      col_idx
    )
  )


  pal <- palettes[[cond]]


  # Draw low-density points first
  # and high-density points last
  ord <- order(
    d$density
  )


  points(
    d$y[ord],
    d$x[ord],
    pch = 16,
    cex = point_size,
    col = pal[
      col_idx[ord]
    ]
  )
}


# ============================================================
# LEGEND
# ============================================================

# legend(
#   "topright",
#   legend = c(
#     "315 + 345",
#     "315 only",
#     "345 only"
#   ),
#   col = c(
#     "grey60",
#     "green",
#     "magenta"
#   ),
#   pch = 16,
#   bty = "n"
# )
quartz.save(file="./figures/pdf/figure01_panelH_scatterplot.pdf", type="pdf")

# ============================================================
# SAMPLE SIZE AND REPLICATE STATISTICS
# ============================================================

# Number of cells per replicate and condition
replicate_counts <- aggregate(
  x ~ condition + file,
  data = data,
  FUN = length
)

names(replicate_counts)[3] <- "n_cells"

cat("\nCells per replicate:\n")
print(replicate_counts, row.names = FALSE)

# Summary statistics per condition
summary_stats <- do.call(
  rbind,
  lapply(split(replicate_counts, replicate_counts$condition), function(d) {
    data.frame(
      condition = unique(d$condition),
      n_replicates = nrow(d),
      total_cells = sum(d$n_cells),
      mean_cells = mean(d$n_cells),
      sd_cells = sd(d$n_cells),
      median_cells = median(d$n_cells),
      min_cells = min(d$n_cells),
      max_cells = max(d$n_cells)
    )
  })
)

rownames(summary_stats) <- NULL

cat("\nSummary per experimental condition:\n")
print(summary_stats, digits = 4, row.names = FALSE)

# Overall summary
cat("\nTotal number of cells:", nrow(data), "\n")
cat("Total number of replicates:", nrow(replicate_counts), "\n")
cat("Mean cells per replicate:", mean(replicate_counts$n_cells), "\n")
cat("SD cells per replicate:", sd(replicate_counts$n_cells), "\n")

# ============================================================
# GENERATE FIGURE LEGEND TEXT
# ============================================================

cat("\nFigure legend statistics:\n")

for (i in seq_len(nrow(summary_stats))) {

  s <- summary_stats[i, ]

  cat(sprintf(
    "%s: n = %d cells across %d replicates (%.1f +/- %.1f cells per replicate, mean +/- SD).\n",
    s$condition,
    s$total_cells,
    s$n_replicates,
    s$mean_cells,
    s$sd_cells
  ))
}
