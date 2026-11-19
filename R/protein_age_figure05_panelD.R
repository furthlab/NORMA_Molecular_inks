# ============================================================
# NORMA molecular inks
# Figure 4D - Nascent-to-total protein ratio
# Nature Methods reporting summary and figure legend
# ============================================================

# ============================================================
# 1. SETTINGS
# ============================================================

in_file <- "./data/csv/protein_age/intensity_per_cell.csv"
out_file <- "./figures/pdf/violin_plot.pdf"

show_plot <- TRUE
save_pdf <- TRUE

fig_width <- 227.0979 / 30
fig_height <- 149.6067 / 30

nascent_ch <- "LAHA_NORMA315"
total_ch <- "NHS_NORMA345"

comp_ids <- c(
  "cytosol",
  "ER",
  "nucleoplasm",
  "nucleolus"
)

comp_labels <- c(
  "Cytosol\n(peripheral)",
  "ER\n(perinuclear)",
  "Nucleoplasm",
  "Nucleoli"
)

tp_palette <- c(
  "#8DB596",
  "#2F5E3F",
  "#6A994E",
  "#7D5BA6"
)

cyt_group <- c("cytosol", "ER")
nuc_group <- c("nucleoplasm", "nucleolus")

min_total <- 1

exclude_fov <- character(0)

plot_ylim <- c(0.45, 2.2)

point_alpha <- 0.25
point_cex <- 0.35
violin_hw <- 0.16

# Preserve the statistical tests from the original script.
# Nucleoli: one-sided Welch t-test.
# Cytoplasm vs nucleus: two-sided paired t-test.

nucleolar_alternative <- "less"

# ============================================================
# 2. LOAD AND VALIDATE DATA
# ============================================================

d <- read.csv(
  in_file,
  stringsAsFactors = FALSE
)

required_columns <- c(
  "timepoint",
  "hours",
  "fov",
  "cell_id",
  "object_type",
  "channel",
  "intensity",
  "n_px"
)

missing_columns <- setdiff(
  required_columns,
  names(d)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

audit <- data.frame(
  Stage = character(),
  Unit = character(),
  Count = integer()
)

record <- function(stage, unit, count) {

  audit <<- rbind(
    audit,
    data.frame(
      Stage = stage,
      Unit = unit,
      Count = count
    )
  )
}

record(
  "Original input",
  "Channel measurements",
  nrow(d)
)

# Keep selected compartments and fluorescence channels

d <- d[
  d$object_type %in% comp_ids &
    d$channel %in% c(nascent_ch, total_ch),
]

# Optional exclusion of specific FOVs

d <- d[
  !paste(d$timepoint, d$fov, sep = "_") %in%
    exclude_fov,
]

record(
  "Selected channels and compartments",
  "Channel measurements",
  nrow(d)
)

# ============================================================
# 3. MATCH NASCENT AND TOTAL CHANNELS
# ============================================================

key <- c(
  "timepoint",
  "hours",
  "fov",
  "cell_id",
  "object_type"
)

if (anyNA(d[, key])) {
  stop("Missing cell or FOV identifiers.")
}

if (anyDuplicated(d[, c(key, "channel")])) {
  stop("Duplicate cell-compartment-channel identifiers.")
}

n <- d[
  d$channel == nascent_ch,
  c(key, "intensity")
]

names(n)[6] <- "nascent"

t <- d[
  d$channel == total_ch,
  c(key, "intensity", "n_px")
]

names(t)[6] <- "total"

r <- merge(n, t, by = key)

record(
  "Matched channel measurements",
  "Cell-compartment pairs",
  nrow(r)
)

record(
  "Unmatched nascent measurements",
  "Channel measurements",
  nrow(n) - nrow(r)
)

record(
  "Unmatched total measurements",
  "Channel measurements",
  nrow(t) - nrow(r)
)

# ============================================================
# 4. QUALITY CONTROL
# ============================================================

valid <- (
  is.finite(r$nascent) &
    is.finite(r$total) &
    is.finite(r$n_px) &
    r$n_px > 0
)

record(
  "Excluded invalid measurements",
  "Cell-compartment pairs",
  sum(!valid)
)

r <- r[valid, ]

n_low <- sum(r$total <= min_total)

record(
  "Excluded total intensity <= 1",
  "Cell-compartment pairs",
  n_low
)

r <- r[r$total > min_total, ]

r$cell_key <- paste(
  r$timepoint,
  r$fov,
  r$cell_id,
  sep = "_"
)

# ============================================================
# 5. WHOLE-CELL NORMALIZATION
# ============================================================

# Calculate pixel-weighted whole-cell fluorescence
# independently for each channel.

rr <- r

rr$wn <- rr$nascent * rr$n_px
rr$wt <- rr$total * rr$n_px

ref <- aggregate(
  cbind(wn, wt, n_px) ~ cell_key + timepoint + fov,
  data = rr,
  FUN = sum
)

ref$ref_nascent <- ref$wn / ref$n_px
ref$ref_total <- ref$wt / ref$n_px

ref <- ref[
  is.finite(ref$ref_nascent) &
    is.finite(ref$ref_total) &
    ref$ref_nascent > 0 &
    ref$ref_total > min_total,
]

n_noref <- length(
  setdiff(
    unique(r$cell_key),
    ref$cell_key
  )
)

record(
  "Excluded cells without usable reference",
  "Unique cells",
  n_noref
)

r <- merge(
  r,
  ref[, c(
    "cell_key",
    "ref_nascent",
    "ref_total"
  )],
  by = "cell_key"
)

r$nascent_norm <- r$nascent / r$ref_nascent
r$total_norm <- r$total / r$ref_total

r$ratio <- r$nascent_norm / r$total_norm

record(
  "Excluded non-finite ratios",
  "Cell-compartment pairs",
  sum(!is.finite(r$ratio))
)

r <- r[is.finite(r$ratio), ]

if (nrow(r) == 0) {
  stop("No observations remain after filtering.")
}

record(
  "Final analysis",
  "Cell-compartment measurements",
  nrow(r)
)

record(
  "Final analysis",
  "Unique cells",
  length(unique(r$cell_key))
)

# ============================================================
# 6. TIMEPOINTS
# ============================================================

tp <- unique(
  r[, c("timepoint", "hours")]
)

tp <- tp[order(tp$hours), ]

tp_levels <- tp$timepoint

tp_cols <- setNames(
  tp_palette[seq_along(tp_levels)],
  tp_levels
)

# ============================================================
# 7. UNIQUE CELL AND FOV COUNTS
# ============================================================

cells <- unique(
  r[, c(
    "timepoint",
    "hours",
    "fov",
    "cell_id",
    "cell_key"
  )]
)

cells_per_fov <- aggregate(
  cell_key ~ timepoint + hours + fov,
  data = cells,
  FUN = length
)

names(cells_per_fov)[4] <- "n_cells"

sample_summary <- do.call(
  rbind,
  lapply(tp_levels, function(tt) {

    z <- cells_per_fov[
      cells_per_fov$timepoint == tt,
    ]

    data.frame(
      timepoint = tt,
      hours = unique(z$hours),
      n_cells = sum(z$n_cells),
      n_fov = nrow(z),
      mean_cells_per_fov = mean(z$n_cells),
      sd_cells_per_fov = sd(z$n_cells),
      median_cells_per_fov = median(z$n_cells),
      min_cells_per_fov = min(z$n_cells),
      max_cells_per_fov = max(z$n_cells)
    )
  })
)

rownames(sample_summary) <- NULL

# ============================================================
# 8. COMPARTMENT SUMMARY STATISTICS
# ============================================================

groups <- split(
  r,
  list(r$timepoint, r$object_type),
  drop = TRUE
)

comp_summary <- do.call(
  rbind,
  lapply(groups, function(z) {

    fm <- tapply(
      z$ratio,
      z$fov,
      median
    )

    data.frame(
      timepoint = z$timepoint[1],
      hours = z$hours[1],
      compartment = z$object_type[1],
      n_cells = length(unique(z$cell_key)),
      n_fov = length(fm),
      mean_cells_per_fov =
        nrow(z) / length(fm),
      median = median(z$ratio),
      q25 = unname(quantile(z$ratio, 0.25)),
      q75 = unname(quantile(z$ratio, 0.75)),
      mean = mean(z$ratio),
      sd = sd(z$ratio),
      fov_median_mean = mean(fm),
      fov_median_sd = sd(fm),
      fov_median_min = min(fm),
      fov_median_max = max(fm)
    )
  })
)

comp_summary <- comp_summary[
  order(
    comp_summary$hours,
    match(comp_summary$compartment, comp_ids)
  ),
]

rownames(comp_summary) <- NULL

# ============================================================
# 9. FOV-LEVEL MEDIANS
# ============================================================

fov_medians <- aggregate(
  ratio ~ timepoint + hours + fov + object_type,
  data = r,
  FUN = median
)

names(fov_medians)[5] <- "median_ratio"

# ============================================================
# 10. NUCLEOLAR WELCH T-TEST
# ============================================================

welch_result <- NULL
welch_test <- NULL

if (length(tp_levels) >= 2) {

  tp1 <- tp_levels[1]
  tp2 <- tp_levels[2]

  a <- fov_medians$median_ratio[
    fov_medians$object_type == "nucleolus" &
      fov_medians$timepoint == tp1
  ]

  b <- fov_medians$median_ratio[
    fov_medians$object_type == "nucleolus" &
      fov_medians$timepoint == tp2
  ]

  if (length(a) >= 2 && length(b) >= 2) {

    welch_test <- tryCatch(
      t.test(
        a,
        b,
        paired = FALSE,
        var.equal = FALSE,
        alternative = nucleolar_alternative,
        conf.level = 0.95
      ),
      error = function(e) {
        message("Welch test unavailable: ", e$message)
        NULL
      }
    )

    if (!is.null(welch_test)) {

      welch_result <- data.frame(
        comparison = "Nucleoli",
        group1 = tp1,
        group2 = tp2,
        n_fov_1 = length(a),
        n_fov_2 = length(b),
        mean1 = mean(a),
        mean2 = mean(b),
        sd1 = sd(a),
        sd2 = sd(b),
        difference_1_minus_2 = mean(a) - mean(b),
        t = unname(welch_test$statistic),
        df = unname(welch_test$parameter),
        p = welch_test$p.value,
        CI_lower = welch_test$conf.int[1],
        CI_upper = welch_test$conf.int[2],
        alternative = nucleolar_alternative,
        correction = "None"
      )
    }
  }
}

# ============================================================
# 11. CYTOPLASM VS NUCLEUS
# ============================================================

# Pool compartments using pixel-weighted channel intensities
# before calculating the regional ratio.

region_ratio <- function(compartments) {

  z <- r[r$object_type %in% compartments, ]

  z$wn <- z$nascent_norm * z$n_px
  z$wt <- z$total_norm * z$n_px

  ag <- aggregate(
    cbind(wn, wt, n_px) ~
      cell_key + timepoint + fov,
    data = z,
    FUN = sum
  )

  ag$ratio <- ag$wn / ag$wt

  ag[, c(
    "cell_key",
    "timepoint",
    "fov",
    "ratio"
  )]
}

pc <- merge(
  region_ratio(cyt_group),
  region_ratio(nuc_group),
  by = c(
    "cell_key",
    "timepoint",
    "fov"
  ),
  suffixes = c("_cyt", "_nuc")
)

pc <- pc[
  is.finite(pc$ratio_cyt) &
    is.finite(pc$ratio_nuc),
]

# Calculate paired FOV medians.

paired_fov <- aggregate(
  cbind(ratio_cyt, ratio_nuc) ~ timepoint + fov,
  data = pc,
  FUN = median
)

paired_fov$difference <-
  paired_fov$ratio_cyt -
  paired_fov$ratio_nuc

paired_result <- NULL
paired_test <- NULL

if (nrow(paired_fov) >= 2) {

  paired_test <- tryCatch(
    t.test(
      paired_fov$ratio_cyt,
      paired_fov$ratio_nuc,
      paired = TRUE,
      alternative = "two.sided",
      conf.level = 0.95
    ),
    error = function(e) {
      message("Paired test unavailable: ", e$message)
      NULL
    }
  )

  if (!is.null(paired_test)) {

    paired_result <- data.frame(
      comparison = "Cytoplasm vs nucleus",
      timepoints = paste(tp_levels, collapse = " + "),
      n_paired_fov = nrow(paired_fov),
      n_paired_cells = nrow(pc),
      mean_cytoplasm = mean(paired_fov$ratio_cyt),
      mean_nucleus = mean(paired_fov$ratio_nuc),
      mean_difference = mean(paired_fov$difference),
      sd_difference = sd(paired_fov$difference),
      t = unname(paired_test$statistic),
      df = unname(paired_test$parameter),
      p = paired_test$p.value,
      CI_lower = paired_test$conf.int[1],
      CI_upper = paired_test$conf.int[2],
      alternative = "two.sided",
      correction = "None"
    )
  }
}

# ============================================================
# 12. NORMALITY ASSESSMENT
# ============================================================

normality_test <- function(label, x) {

  x <- x[is.finite(x)]

  if (
    length(x) < 3 ||
    length(x) > 5000 ||
    length(unique(x)) < 2
  ) {

    return(
      data.frame(
        group = label,
        n = length(x),
        W = NA_real_,
        p = NA_real_,
        note = "Not evaluable"
      )
    )
  }

  st <- shapiro.test(x)

  data.frame(
    group = label,
    n = length(x),
    W = unname(st$statistic),
    p = st$p.value,
    note = if (length(x) < 10) {
      "Limited power at small sample size"
    } else {
      ""
    }
  )
}

normality <- do.call(
  rbind,
  lapply(tp_levels, function(tt) {

    normality_test(
      paste("Nucleoli", tt),
      fov_medians$median_ratio[
        fov_medians$object_type == "nucleolus" &
          fov_medians$timepoint == tt
      ]
    )
  })
)

normality <- rbind(
  normality,
  normality_test(
    "Paired cytoplasm-nucleus differences",
    paired_fov$difference
  )
)

# ============================================================
# 13. OFF-SCALE OBSERVATIONS
# ============================================================

# These observations are not excluded from statistical analysis.

offscale <- do.call(
  rbind,
  lapply(groups, function(z) {

    data.frame(
      timepoint = z$timepoint[1],
      compartment = z$object_type[1],
      below = sum(z$ratio < plot_ylim[1]),
      above = sum(z$ratio > plot_ylim[2])
    )
  })
)

rownames(offscale) <- NULL

# ============================================================
# 14. VIOLIN PLOT FUNCTIONS
# ============================================================

violin_shape <- function(v, hw) {

  dd <- density(v, n = 512)

  inside <- (
    dd$x >= min(v) &
      dd$x <= max(v)
  )

  list(
    y = dd$x[inside],
    w = dd$y[inside] / max(dd$y[inside]) * hw
  )
}

draw_violin <- function(shape, x, col) {

  polygon(
    c(x - shape$w, rev(x + shape$w)),
    c(shape$y, rev(shape$y)),
    col = col,
    border = NA
  )
}

sina_x <- function(v, shape, x) {

  w <- approx(
    shape$y,
    shape$w,
    xout = v,
    rule = 2
  )$y

  x + runif(length(v), -1, 1) * w * 0.95
}

format_p <- function(p) {

  if (p < 0.001) {
    formatC(p, format = "e", digits = 2)
  } else {
    formatC(p, format = "f", digits = 4)
  }
}

# ============================================================
# 15. DRAW FIGURE
# ============================================================

draw_figure <- function() {

  nt <- length(tp_levels)

  offsets <- (
    seq_len(nt) - (nt + 1) / 2
  ) * 0.36

  ylim <- plot_ylim

  par(
    mar = c(4.5, 5.5, 3, 1),
    mgp = c(3.4, 0.7, 0),
    oma = c(3.4, 0, 0, 0)
  )

  plot(
    NA,
    xlim = c(0.5, length(comp_ids) + 0.5),
    ylim = ylim,
    axes = FALSE,
    xlab = "",
    ylab = paste0(
      "Nascent / total, each channel\n",
      "normalized to the cell's own whole cell"
    )
  )

  abline(
    h = pretty(ylim),
    col = "grey92",
    lwd = 0.6
  )

  abline(
    h = 1,
    lty = 2,
    col = "grey40"
  )

  axis(2, las = 1)

  set.seed(1)

  for (k in seq_len(nt)) {

    meds <- rep(NA_real_, length(comp_ids))

    for (j in seq_along(comp_ids)) {

      v <- r$ratio[
        r$object_type == comp_ids[j] &
          r$timepoint == tp_levels[k]
      ]

      if (!length(v)) next

      x <- j + offsets[k]

      hi <- v > ylim[2]
      lo <- v < ylim[1]

      vin <- v[!hi & !lo]

      # Violin and individual cells

      if (
        length(vin) >= 3 &&
        diff(range(vin)) > 0
      ) {

        shape <- violin_shape(vin, violin_hw)

        draw_violin(
          shape,
          x,
          adjustcolor(tp_cols[k], 0.18)
        )

        points(
          sina_x(vin, shape, x),
          vin,
          pch = 16,
          cex = point_cex,
          col = adjustcolor(
            tp_cols[k],
            point_alpha
          )
        )

      } else if (length(vin)) {

        points(
          rep(x, length(vin)),
          vin,
          pch = 16,
          cex = point_cex,
          col = tp_cols[k]
        )
      }

      # Off-scale observations

      if (any(hi)) {

        points(
          x + runif(sum(hi), -0.05, 0.05),
          rep(ylim[2], sum(hi)),
          pch = 17,
          cex = 0.5,
          col = tp_cols[k],
          xpd = NA
        )
      }

      if (any(lo)) {

        points(
          x + runif(sum(lo), -0.05, 0.05),
          rep(ylim[1], sum(lo)),
          pch = 25,
          cex = 0.5,
          col = tp_cols[k],
          bg = tp_cols[k],
          xpd = NA
        )
      }

      # Boxplot includes all observations.

      boxplot(
        v,
        at = x,
        add = TRUE,
        boxwex = 0.07,
        outline = FALSE,
        axes = FALSE,
        col = adjustcolor("white", 0.8),
        lwd = 0.9,
        staplewex = 0,
        whisklty = 1
      )

      text(
        x,
        ylim[1] + 0.02 * diff(ylim),
        length(v),
        cex = 0.6,
        col = "grey35",
        adj = c(0.5, 0)
      )

      meds[j] <- median(v)
    }

    lines(
      seq_along(comp_ids) + offsets[k],
      meds,
      col = tp_cols[k],
      lwd = 1.6
    )
  }

  # ----------------------------------------------------------
  # Nucleolar significance bracket
  # ----------------------------------------------------------

  if (!is.null(welch_result)) {

    j <- match("nucleolus", comp_ids)

    x1 <- j + offsets[1]
    x2 <- j + offsets[2]

    v <- r$ratio[
      r$object_type == "nucleolus" &
        r$timepoint %in% tp_levels[1:2]
    ]

    v <- v[v <= ylim[2]]

    if (length(v)) {

      yb <- min(
        max(v) + 0.05 * diff(ylim),
        ylim[2] - 0.10 * diff(ylim)
      )

      tk <- 0.015 * diff(ylim)

      segments(
        c(x1, x1, x2),
        c(yb - tk, yb, yb),
        c(x1, x2, x2),
        c(yb, yb, yb - tk),
        lwd = 0.9
      )

      text(
        mean(c(x1, x2)),
        yb + 0.012 * diff(ylim),
        paste0(
          "p = ", format_p(welch_result$p),
          "\n(n = ", welch_result$n_fov_1,
          " vs ", welch_result$n_fov_2,
          " FOVs)"
        ),
        adj = c(0.5, 0),
        cex = 0.7
      )
    }
  }

  # ----------------------------------------------------------
  # Cytoplasm vs nucleus significance bracket
  # ----------------------------------------------------------

  if (!is.null(paired_result)) {

    xc <- range(match(cyt_group, comp_ids)) +
      c(-0.25, 0.25)

    xn <- range(match(nuc_group, comp_ids)) +
      c(-0.25, 0.25)

    v <- r$ratio[r$ratio <= ylim[2]]

    if (length(v)) {

      yb <- min(
        unname(quantile(v, 0.995)) +
          0.03 * diff(ylim),
        ylim[2] - 0.08 * diff(ylim)
      )

      tk <- 0.012 * diff(ylim)
      rise <- 0.03 * diff(ylim)

      cx <- mean(xc)
      nx <- mean(xn)

      segments(
        c(xc[1], xn[1]),
        yb,
        c(xc[2], xn[2]),
        yb,
        col = "grey20",
        lwd = 1.1
      )

      segments(
        c(xc, xn),
        yb,
        c(xc, xn),
        yb - tk,
        col = "grey20",
        lwd = 1.1
      )

      segments(
        c(cx, nx),
        yb,
        c(cx, nx),
        yb + rise,
        col = "grey20",
        lwd = 1.1
      )

      segments(
        cx,
        yb + rise,
        nx,
        yb + rise,
        col = "grey20",
        lwd = 1.1
      )

      lab <- paste0(
        "p = ", format_p(paired_result$p),
        " (paired, n = ",
        paired_result$n_paired_fov,
        " FOVs)"
      )

      xm <- mean(c(cx, nx))
      yt <- yb + rise

      w <- strwidth(lab, cex = 0.68) * 1.08
      h <- strheight(lab, cex = 0.68) * 1.6

      rect(
        xm - w / 2,
        yt - h / 2,
        xm + w / 2,
        yt + h / 2,
        col = "white",
        border = NA
      )

      text(
        xm,
        yt,
        lab,
        cex = 0.68,
        col = "grey20"
      )
    }
  }

  # ----------------------------------------------------------
  # Axes and legend
  # ----------------------------------------------------------

  axis(
    1,
    at = seq_along(comp_ids),
    labels = FALSE
  )

  mtext(
    comp_labels,
    side = 1,
    at = seq_along(comp_ids),
    line = 2,
    cex = 0.8
  )

  box(bty = "l")

  legend(
    "topright",
    legend = tp_levels,
    pch = 15,
    pt.cex = 1.5,
    col = adjustcolor(tp_cols, 0.6),
    bty = "n",
    horiz = TRUE,
    cex = 0.9
  )

  mtext(
    paste0(
      "Nascent-to-total ratio relative to the whole cell",
      " (dashed line = whole cell)"
    ),
    side = 3,
    line = 1,
    adj = 0,
    font = 2,
    cex = 1.05
  )

  mtext(
    paste0(
      "Violin = distribution of cells; ",
      "points = individual cells; ",
      "boxes = quartiles; ",
      "lines join medians; numbers = n cells"
    ),
    side = 1,
    outer = TRUE,
    adj = 0.02,
    cex = 0.75,
    col = "grey30",
    line = 0.2
  )

  n_out <- sum(
    r$ratio < ylim[1] |
      r$ratio > ylim[2]
  )

  mtext(
    paste(
      n_out,
      "cell-compartment measurements off-scale (triangles)"
    ),
    side = 1,
    outer = TRUE,
    adj = 0.02,
    cex = 0.75,
    col = "grey30",
    line = 1.2
  )

  mtext(
    paste0(
      "Nucleoli: ",
      if (nucleolar_alternative == "two.sided") {
        "two-sided"
      } else {
        "one-sided"
      },
      " Welch t-test; cytoplasm vs nucleus: ",
      "two-sided paired t-test; both on FOV medians"
    ),
    side = 1,
    outer = TRUE,
    adj = 0.02,
    cex = 0.75,
    col = "grey30",
    line = 2.2
  )
}

# ============================================================
# 16. QUARTZ WINDOW
# ============================================================

if (show_plot && capabilities("aqua")) {

  quartz(
    title = "Figure 4D",
    width = fig_width,
    height = fig_height
  )

  draw_figure()
}

# ============================================================
# 17. SAVE MANUSCRIPT FIGURE
# ============================================================

if (save_pdf) {

  dir.create(
    dirname(out_file),
    recursive = TRUE,
    showWarnings = FALSE
  )

  pdf(
    out_file,
    width = fig_width,
    height = fig_height
  )

  draw_figure()

  invisible(dev.off())
}

# ============================================================
# 18. NATURE METHODS CONSOLE REPORT
# ============================================================

old_options <- options(digits = 16)

cat("\n")
cat("============================================================\n")
cat("FIGURE 4D - NATURE METHODS REPORTING SUMMARY\n")
cat("============================================================\n")

cat("\nR version:\n")
cat(R.version.string, "\n")

cat("\nTotal unique cells:", nrow(cells), "\n")
cat("Total cell-compartment measurements:", nrow(r), "\n")
cat("Total FOVs:", nrow(cells_per_fov), "\n")

cat("\nOverall number of cells per FOV:\n")

cat(
  sprintf(
    "%.2f +/- %.2f cells (mean +/- SD)\n",
    mean(cells_per_fov$n_cells),
    sd(cells_per_fov$n_cells)
  )
)

cat("\n------------------------------------------------------------\n")
cat("SAMPLE SIZES BY TIMEPOINT\n")
cat("------------------------------------------------------------\n\n")

print(
  sample_summary,
  row.names = FALSE
)

cat("\nSample size descriptions:\n")

for (i in seq_len(nrow(sample_summary))) {

  z <- sample_summary[i, ]

  cat(
    sprintf(
      "%g h: n = %d cells across %d FOVs (%.1f +/- %.1f cells per FOV, mean +/- SD).\n",
      z$hours,
      z$n_cells,
      z$n_fov,
      z$mean_cells_per_fov,
      z$sd_cells_per_fov
    )
  )
}

cat("\n------------------------------------------------------------\n")
cat("NUMBER OF CELLS PER FOV\n")
cat("------------------------------------------------------------\n\n")

print(
  cells_per_fov,
  row.names = FALSE
)

cat("\n------------------------------------------------------------\n")
cat("COMPARTMENT SUMMARY\n")
cat("------------------------------------------------------------\n\n")

print(
  comp_summary,
  row.names = FALSE
)

cat("\n------------------------------------------------------------\n")
cat("NUCLEOLAR WELCH T-TEST\n")
cat("------------------------------------------------------------\n\n")

cat("Experimental unit: FOV median\n")
cat("Alternative:", nucleolar_alternative, "\n")
cat("Multiple-comparison correction: None\n\n")

if (!is.null(welch_result)) {

  print(
    welch_result,
    row.names = FALSE
  )

  print(welch_test)

} else {

  cat("Insufficient data or test unavailable.\n")
}

cat("\n------------------------------------------------------------\n")
cat("CYTOPLASM VS NUCLEUS PAIRED T-TEST\n")
cat("------------------------------------------------------------\n\n")

cat("Experimental unit: Paired FOV medians\n")
cat("Test: Two-sided paired t-test\n")
cat("Timepoints pooled\n")
cat("Multiple-comparison correction: None\n\n")

if (!is.null(paired_result)) {

  print(
    paired_result,
    row.names = FALSE
  )

  print(paired_test)

} else {

  cat("Insufficient data or test unavailable.\n")
}

cat("\nPaired FOV medians:\n")

print(
  paired_fov,
  row.names = FALSE
)

cat("\n------------------------------------------------------------\n")
cat("SHAPIRO-WILK NORMALITY TESTS\n")
cat("------------------------------------------------------------\n\n")

print(
  normality,
  row.names = FALSE
)

cat(
  "\nNormality tests were performed on FOV-level values",
  "and paired differences, as appropriate.\n"
)

cat(
  "Interpretation is limited by the small number of FOVs.\n"
)

cat("\n------------------------------------------------------------\n")
cat("DATA EXCLUSION AUDIT\n")
cat("------------------------------------------------------------\n\n")

print(
  audit,
  row.names = FALSE
)

cat("\nExplicitly excluded FOVs: ")

if (length(exclude_fov)) {

  cat(
    paste(exclude_fov, collapse = ", "),
    "\n"
  )

} else {

  cat("None\n")
}

cat("\n------------------------------------------------------------\n")
cat("OFF-SCALE OBSERVATIONS\n")
cat("------------------------------------------------------------\n\n")

print(
  offscale,
  row.names = FALSE
)

cat(
  "\nOff-scale observations were retained in all",
  "statistical analyses.\n"
)

cat("\n------------------------------------------------------------\n")
cat("EXPERIMENTAL DESIGN\n")
cat("------------------------------------------------------------\n\n")

cat("Individual cells are observational units.\n")
cat("Statistical tests are based on FOV medians.\n")
cat("Biological replicate counts cannot be inferred from FOV IDs.\n")
cat("No multiple-comparison correction was performed.\n")

cat("\n------------------------------------------------------------\n")
cat("SOFTWARE\n")
cat("------------------------------------------------------------\n\n")

print(sessionInfo())

options(old_options)

# ============================================================
# 19. AUTOMATIC LATEX FIGURE LEGEND
# ============================================================

latex_p <- function(p) {

  if (!is.finite(p)) {
    return("P value unavailable")
  }

  if (p == 0) {
    return("$P$ below numerical precision")
  }

  if (p < 0.001) {

    z <- strsplit(
      sprintf("%.3e", p),
      "e",
      fixed = TRUE
    )[[1]]

    return(
      sprintf(
        "$P = %s \\times 10^{%d}$",
        z[1],
        as.integer(z[2])
      )
    )
  }

  sprintf("$P = %.4g$", p)
}

# ------------------------------------------------------------
# Sample size description
# ------------------------------------------------------------

cell_text <- vapply(
  seq_len(nrow(sample_summary)),
  function(i) {

    z <- sample_summary[i, ]

    sd_text <- if (is.finite(z$sd_cells_per_fov)) {

      sprintf(
        "%.1f $\\pm$ %.1f cells per FOV, mean $\\pm$ SD",
        z$mean_cells_per_fov,
        z$sd_cells_per_fov
      )

    } else {

      sprintf(
        "%d cells in one FOV",
        z$n_cells
      )
    }

    sprintf(
      "%d cells across %d FOVs after a %g~h pulse (%s)",
      z$n_cells,
      z$n_fov,
      z$hours,
      sd_text
    )
  },
  character(1)
)

# ------------------------------------------------------------
# Compartment-specific cell counts
# ------------------------------------------------------------

comp_text <- vapply(
  seq_len(nrow(sample_summary)),
  function(i) {

    z <- sample_summary[i, ]

    cc <- comp_summary[
      comp_summary$timepoint == z$timepoint,
    ]

    counts <- cc$n_cells[
      match(comp_ids, cc$compartment)
    ]

    counts[is.na(counts)] <- 0

    sprintf(
      "%g~h: %s",
      z$hours,
      paste(counts, collapse = ", ")
    )
  },
  character(1)
)

# ------------------------------------------------------------
# Statistical test description
# ------------------------------------------------------------

stat_text <- character(0)

if (!is.null(welch_result)) {

  h1 <- sample_summary$hours[
    match(
      welch_result$group1,
      sample_summary$timepoint
    )
  ]

  h2 <- sample_summary$hours[
    match(
      welch_result$group2,
      sample_summary$timepoint
    )
  ]

  sidedness <- if (
    nucleolar_alternative == "two.sided"
  ) {
    "two-sided"
  } else {
    "one-sided"
  }

  direction <- switch(
    nucleolar_alternative,
    less = sprintf(
      ", testing whether the ratio at %g~h was lower than at %g~h",
      h1, h2
    ),
    greater = sprintf(
      ", testing whether the ratio at %g~h was higher than at %g~h",
      h1, h2
    ),
    two.sided = ""
  )

  stat_text <- c(
    stat_text,
    sprintf(
      paste0(
        "Nucleolar ratios at %g and %g~h ",
        "were compared using a %s Welch's ",
        "$t$-test on FOV medians%s ",
        "($n = %d$ and $%d$ FOVs, respectively; %s)."
      ),
      h1,
      h2,
      sidedness,
      direction,
      welch_result$n_fov_1,
      welch_result$n_fov_2,
      latex_p(welch_result$p)
    )
  )
}

if (!is.null(paired_result)) {

  stat_text <- c(
    stat_text,
    sprintf(
      paste0(
        "Cytoplasmic (cytosol and perinuclear region) ",
        "and nuclear (nucleoplasm and nucleoli) ratios ",
        "were compared using a two-sided paired ",
        "$t$-test on FOV medians, pooling pulse durations ",
        "($n = %d$ paired FOVs; %s)."
      ),
      paired_result$n_paired_fov,
      latex_p(paired_result$p)
    )
  )
}

# ------------------------------------------------------------
# Construct panel D legend
# ------------------------------------------------------------

panel_d <- paste(

  "\\textbf{d)} Nascent-to-total protein ratio per cell and compartment, calculated after normalizing each fluorescence channel to the pixel-weighted whole-cell mean (dashed line, ratio of 1).",

  paste0(
    "Data comprise ",
    paste(cell_text, collapse = "; "),
    "."
  ),

  paste0(
    "Cell counts for cytosol, perinuclear region, ",
    "nucleoplasm and nucleoli, respectively: ",
    paste(comp_text, collapse = "; "),
    "."
  ),

  paste0(
    "Points represent individual cell-compartment ",
    "measurements; violin widths indicate kernel ",
    "density estimates of values within the ",
    "displayed range. Boxes indicate the median ",
    "and approximate interquartile range, with ",
    "whiskers extending to the most extreme ",
    "observations within 1.5 times the ",
    "interquartile range. Connecting lines join ",
    "compartment-specific medians."
  ),

  if (sum(offscale$below + offscale$above) > 0) {

    paste0(
      "Triangles indicate off-scale observations ",
      "retained in the statistical analysis."
    )

  } else {
    ""
  },

  paste(stat_text, collapse = " "),

  if (length(stat_text)) {
    "P values were not adjusted for multiple comparisons."
  } else {
    ""
  }
)

# ------------------------------------------------------------
# Complete figure legend
# ------------------------------------------------------------

figure_legend <- c(

  "\\begin{figure*}[t]",

  "    \\centering",

  paste0(
    "    \\includegraphics[width=\\textwidth,",
    "height=0.60\\textheight,keepaspectratio]",
    "{figures/NORMA_Nature_Methods_Figure04.pdf}"
  ),

  "\\caption{",

  paste0(
    "\\textbf{\\textcolor{blue}{",
    "Molecular inks map protein age within single cells.}}"
  ),

  paste0(
    "\\textbf{a)} Metabolic labeling strategy. ",
    "L-Azidohomoalanine (L-AHA) was incorporated ",
    "into nascent proteins and ligated to NORMA-315 ",
    "by SPAAC, installing a latent photocaged ",
    "cyclooctyne for subsequent photoactivation ",
    "and fluorophore conjugation."
  ),

  paste0(
    "\\textbf{b)} Labeling scheme: NHS-NORMA-345 ",
    "functionalizes all protein present at fixation, ",
    "whereas NORMA-315 is ligated to proteins ",
    "synthesized during a 2 or 24~h L-AHA pulse."
  ),

  paste0(
    "\\textbf{c)} In-gel fluorescence of metabolically ",
    "labeled proteins. Lane 1, fluorophore only; ",
    "lane 2, NORMA-315 without UV irradiation; ",
    "lane 3, NORMA-315 following 308 nm photoactivation; ",
    "lane 4, DBCO-MB660 positive control. ",
    "Sequential photo-SPAAC labels a broadly similar ",
    "nascent proteome, with increased band smearing ",
    "following prolonged UV irradiation."
  ),

  panel_d,

  paste0(
    "\\textbf{e--g)} U2OS cells: merge, ",
    "total protein (magenta), nascent protein ",
    "(green), DAPI (blue). Arrowheads, e, ",
    "perinuclear zone; f,g, mitotic cells. ",
    "Scale bars, 40~\\textmu m (e), ",
    "20~\\textmu m (f,g)."
  ),

  "}",

  "    \\label{fgr:fig04}",

  "\\end{figure*}"
)

cat("\n\n")
cat("============================================================\n")
cat("COPY-READY LATEX FIGURE 4 LEGEND\n")
cat("============================================================\n\n")

cat(
  paste(figure_legend, collapse = "\n"),
  "\n"
)

cat("\n============================================================\n")
cat("END LATEX FIGURE LEGEND\n")
cat("============================================================\n")

# ============================================================
# 20. REMAINING REPORTING REQUIREMENTS
# ============================================================

cat("\nREMAINING ITEMS TO VERIFY MANUALLY:\n\n")

cat(
  "[ ] Number of independently prepared biological replicates.\n"
)

cat(
  "[ ] Independence of FOVs used for statistical inference.\n"
)

cat(
  "[ ] Number of independent experiments for panel C.\n"
)

cat(
  "[ ] Number of independent experiments for panels E-G.\n"
)

if (nucleolar_alternative != "two.sided") {

  cat(
    "[ ] Justification of the one-sided nucleolar Welch t-test.\n"
  )

  cat(
    "[ ] Update Methods statement that all tests were two-sided.\n"
  )
}

cat("\nAnalysis complete.\n")

if (save_pdf) {
  cat("Figure saved to:", out_file, "\n")
}


cat("Below 0.45:", sum(r$ratio < 0.45), "\n")
cat("Above 2.2:", sum(r$ratio > 2.2), "\n")

cat("Total off-scale measurements:",
    sum(r$ratio < 0.45 | r$ratio > 2.2), "\n")

cat("Unique cells with off-scale measurements:",
    length(unique(r$cell_key[r$ratio < 0.45 | r$ratio > 2.2])), "\n")

cat("Percentage of measurements off-scale:",
    round(100 * mean(r$ratio < 0.45 | r$ratio > 2.2), 2), "%\n")

# Breakdown by pulse duration and compartment
with(r[r$ratio < 0.45 | r$ratio > 2.2, ],
     table(timepoint, object_type))





