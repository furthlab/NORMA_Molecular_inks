# ============================================================
# NORMA spatial proteomics
# Selected subcellular marker proteins
# With UniProt subcellular localization
# ============================================================


# 1. Read TSV files ------------------------------------------

read_sample <- function(file) {

  dat <- read.delim(
    file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  intensity_col <- grep(
    "PG.MS2Quantity$",
    colnames(dat),
    value = TRUE
  )

  stopifnot(length(intensity_col) == 1)

  out <- data.frame(
    Gene = trimws(dat$PG.Genes),
    Intensity = as.numeric(dat[[intensity_col]])
  )

  # Remove ambiguous gene assignments
  out <- subset(
    out,
    !is.na(Gene) &
      Gene != "" &
      !grepl(";", Gene) &
      is.finite(Intensity) &
      Intensity > 0
  )

  # Collapse duplicate gene symbols
  aggregate(Intensity ~ Gene, data = out, FUN = sum)
}

A <- read_sample("./data/csv/LC-MSMS/20261111_Norma_Sample_A.tsv")
B <- read_sample("./data/csv/LC-MSMS/20261111_Norma_Sample_B.tsv")
C <- read_sample("./data/csv/LC-MSMS/20261111_NORMA_Sample_C.tsv")

colnames(A)[2] <- "Control"
colnames(B)[2] <- "Nucleus"
colnames(C)[2] <- "Cytoplasm"


# 2. Merge samples -------------------------------------------

dat <- Reduce(
  function(x, y) merge(x, y, by = "Gene", all = TRUE),
  list(A, B, C)
)


# 3. Calculate median normalization --------------------------

dat$log2FC_raw <- log2(dat$Nucleus / dat$Cytoplasm)

shared <- is.finite(dat$log2FC_raw)

stopifnot(any(shared))

# Calculate normalization before imputation
median_shift <- median(dat$log2FC_raw[shared])

cat("Shared genes:", sum(shared), "\n")
cat("Median shift:", median_shift, "\n")


# 4. Missing value imputation -------------------------------

dat$Nucleus_detected <- is.finite(dat$Nucleus)
dat$Cytoplasm_detected <- is.finite(dat$Cytoplasm)

# Identify proteins with one missing compartment
dat$Imputed <- xor(
  dat$Nucleus_detected,
  dat$Cytoplasm_detected
)

# Imputation parameters
imputation_quantile <- 0.01
imputation_fraction <- 0.5

# Sample-specific low-intensity estimates
floor_B <- imputation_fraction * quantile(
  dat$Nucleus[dat$Nucleus_detected],
  probs = imputation_quantile,
  names = FALSE
)

floor_C <- imputation_fraction * quantile(
  dat$Cytoplasm[dat$Cytoplasm_detected],
  probs = imputation_quantile,
  names = FALSE
)

# Preserve original measurements
dat$Nucleus_imp <- dat$Nucleus
dat$Cytoplasm_imp <- dat$Cytoplasm

# Impute missing nucleus measurements
idx_B <- !dat$Nucleus_detected &
  dat$Cytoplasm_detected

dat$Nucleus_imp[idx_B] <- floor_B

# Impute missing cytoplasm measurements
idx_C <- dat$Nucleus_detected &
  !dat$Cytoplasm_detected

dat$Cytoplasm_imp[idx_C] <- floor_C

# Calculate median-centered log2 fold change
dat$log2FC <- log2(
  dat$Nucleus_imp / dat$Cytoplasm_imp
) - median_shift

# Record detection category
dat$Detection <- ifelse(
  dat$Nucleus_detected & dat$Cytoplasm_detected,
  "Both compartments",
  ifelse(
    dat$Nucleus_detected,
    "Nucleus only",
    ifelse(
      dat$Cytoplasm_detected,
      "Cytoplasm only",
      "Control only"
    )
  )
)

cat("Nucleus imputation:", floor_B, "\n")
cat("Cytoplasm imputation:", floor_C, "\n")

print(table(dat$Detection))


# 5. Select subcellular marker proteins ----------------------

markers <- c(
  "HNRNPC",
  "NPM1",
  "LMNA",
  "LMNB1",
  "SRSF1",
  "HNRNPA1",
  "H4C1",
  "GAPDH",
  "HSPA5",
  "ATP5F1A",
  "PSAP",
  "CTSD",
  "LAMP1",
  "HEXB",
  "TUBB4A"
)

res <- dat[match(markers, dat$Gene), ]

# Remove proteins absent from both compartments
res <- res[is.finite(res$log2FC), ]

# Sort by compartment enrichment
res <- res[order(res$log2FC, decreasing = TRUE), ]

stopifnot(nrow(res) > 0)


# 6. Retrieve subcellular localization from UniProt ----------

fetch_uniprot <- function(genes) {

  genes <- unique(genes[!is.na(genes)])

  gene_query <- paste(
    paste0("gene_exact:", genes),
    collapse = " OR "
  )

  query <- paste0(
    "(organism_id:9606) AND ",
    "(reviewed:true) AND (",
    gene_query,
    ")"
  )

  request_url <- paste0(
    "https://rest.uniprot.org/uniprotkb/search?",
    "query=", URLencode(query, reserved = TRUE),
    "&format=tsv",
    "&fields=accession,gene_primary,gene_names,",
    "cc_subcellular_location",
    "&size=500"
  )

  out <- read.delim(
    request_url,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )

  colnames(out) <- c(
    "UniProt",
    "Gene",
    "Gene_names",
    "Localization_raw"
  )

  return(out)
}

options(timeout = 120)

uniprot <- tryCatch(
  fetch_uniprot(res$Gene),
  error = function(e) {

    warning(
      "UniProt retrieval failed: ",
      conditionMessage(e)
    )

    data.frame(
      UniProt = character(0),
      Gene = character(0),
      Gene_names = character(0),
      Localization_raw = character(0)
    )
  }
)

cat("Retrieved UniProt entries:", nrow(uniprot), "\n")


# 7. Clean and standardize localization annotations ----------

clean_localization <- function(x, max_locations = 2) {

  if (is.na(x) || !nzchar(trimws(x))) {
    return("Not annotated")
  }

  # Remove evidence codes
  x <- gsub("\\{[^}]*\\}", "", x)

  # Remove annotation headings
  x <- gsub(
    "SUBCELLULAR LOCATION:",
    "",
    x,
    fixed = TRUE
  )

  # Remove isoform and protein-product headings
  x <- gsub("\\[[^]]*\\]:", "", x)

  # Remove explanatory notes
  x <- sub("Note=.*$", "", x)

  # Split localization entries
  loc <- unlist(strsplit(x, "\\."))

  # Remove membrane topology descriptions
  loc <- sub(";.*$", "", loc)

  # Remove leading punctuation
  loc <- gsub("^[[:space:];,]+", "", loc)

  # Remove excess whitespace
  loc <- gsub("[[:space:]]+", " ", loc)
  loc <- trimws(loc)

  loc <- loc[nzchar(loc)]

  # Standardize and abbreviate terminology
  replacements <- c(
    "Nucleus lamina" = "Nuclear lamina",
    "Nucleus envelope" = "Nuclear envelope",
    "Nucleus membrane" = "Nuclear membrane",
    "Nucleus, nucleolus" = "Nucleolus",
    "Nucleus, nucleoplasm" = "Nucleoplasm",
    "Nucleus speckle" = "Nuclear speckles",
    "Nucleus, speckle" = "Nuclear speckles",
    "Cytoplasm, cytosol" = "Cytosol",
    "Cytoplasm, cytoskeleton" = "Cytoskeleton",
    "Endoplasmic reticulum lumen" = "ER lumen",
    "Endoplasmic reticulum membrane" = "ER membrane",
    "Endoplasmic reticulum" = "ER",
    "Lysosome membrane" = "Lysosomal membrane",
    "Endosome membrane" = "Endosomal membrane",
    "Mitochondrion inner membrane" = "Mitochondrial inner membrane",
    "Mitochondrion outer membrane" = "Mitochondrial outer membrane",
    "Mitochondrion matrix" = "Mitochondrial matrix",
    "Golgi apparatus membrane" = "Golgi membrane",
    "Golgi apparatus" = "Golgi",
    "Cell membrane" = "Plasma membrane"
  )

  for (term in names(replacements)) {
    loc[loc == term] <- replacements[[term]]
  }

  # Remove duplicate annotations
  loc <- unique(loc)

  if (length(loc) == 0) {
    return("Not annotated")
  }

  # Keep at most two locations
  paste(
    head(loc, max_locations),
    collapse = " / "
  )
}

uniprot$Localization <- vapply(
  uniprot$Localization_raw,
  clean_localization,
  FUN.VALUE = character(1)
)


# 8. Match UniProt entries to marker proteins ----------------

# Match primary gene names first, then alternative names
match_uniprot <- function(gene) {

  idx <- which(uniprot$Gene == gene)

  if (length(idx) == 0) {

    idx <- which(vapply(
      strsplit(
        uniprot$Gene_names,
        "[[:space:]]+"
      ),
      function(x) gene %in% x,
      FUN.VALUE = logical(1)
    ))
  }

  return(idx)
}

# UniProt accession numbers
res$UniProt <- vapply(
  res$Gene,
  function(gene) {

    idx <- match_uniprot(gene)

    if (length(idx) == 0) {
      return(NA_character_)
    }

    paste(
      unique(uniprot$UniProt[idx]),
      collapse = ";"
    )
  },
  FUN.VALUE = character(1)
)

# Preserve complete UniProt annotations
res$Localization_raw <- vapply(
  res$Gene,
  function(gene) {

    idx <- match_uniprot(gene)

    if (length(idx) == 0) {
      return(NA_character_)
    }

    paste(
      unique(uniprot$Localization_raw[idx]),
      collapse = " | "
    )
  },
  FUN.VALUE = character(1)
)

# Generate concise localization annotations
res$Localization <- vapply(
  res$Gene,
  function(gene) {

    idx <- match_uniprot(gene)

    if (length(idx) == 0) {
      return("Not annotated")
    }

    loc <- unique(uniprot$Localization[idx])
    loc <- loc[loc != "Not annotated"]

    if (length(loc) == 0) {
      return("Not annotated")
    }

    loc <- unlist(
      strsplit(loc, " / ", fixed = TRUE)
    )

    paste(
      head(unique(loc), 2),
      collapse = " / "
    )
  },
  FUN.VALUE = character(1)
)

# Print annotated marker proteins
print(
  res[, c(
    "Gene",
    "UniProt",
    "Localization",
    "Control",
    "Nucleus",
    "Cytoplasm",
    "log2FC",
    "Imputed",
    "Detection"
  )],
  row.names = FALSE
)


# 9. Prepare plotting variables -----------------------------

col_nucleus   <- "#5279CB"
col_cytoplasm <- "#DB9251"

plot_cols <- ifelse(
  res$log2FC >= 0,
  col_nucleus,
  col_cytoplasm
)

# Reverse plotting order
values   <- rev(res$log2FC)
genes    <- rev(res$Gene)
location <- rev(res$Localization)
cols     <- rev(plot_cols)
imputed  <- rev(res$Imputed)

# Ensure localization labels occupy only one line
location <- gsub("[\r\n]+", " ", location)

# Limit localization descriptions to two locations
location <- vapply(
  strsplit(location, " / ", fixed = TRUE),
  function(x) {
    paste(
      head(unique(x), 2),
      collapse = " / "
    )
  },
  FUN.VALUE = character(1)
)

# Remove secondary localization if label is too long
max_chars <- 38

location <- vapply(
  location,
  function(x) {

    if (nchar(x) > max_chars) {
      x <- strsplit(x, " / ", fixed = TRUE)[[1]][1]
    }

    if (nchar(x) > max_chars) {
      x <- paste0(
        substr(x, 1, max_chars - 3),
        "..."
      )
    }

    return(x)
  },
  FUN.VALUE = character(1)
)

# Automatically determine symmetric plotting limits
axis_limit <- max(1, ceiling(max(abs(values))))
plot_limit <- axis_limit + 0.9

xlim <- c(-plot_limit, plot_limit)


# 10. Plot using base R --------------------------------------

quartz(width = 11, height = 8)

par(
  family = "sans",
  las = 1,
  bty = "n",
  bg = "white",
  oma = c(0, 0, 0, 0)
)

# More compact annotation panel
# Adjust these values to control horizontal proportions
layout(
  matrix(c(1, 2), nrow = 1),
  widths = c(0.32, 0.68)
)

# Calculate bar positions before plotting
bp <- as.numeric(
  barplot(
    values,
    horiz = TRUE,
    width = 0.7,
    space = 0.65,
    plot = FALSE
  )
)

ylim <- c(
  min(bp) - 0.8,
  max(bp) + 0.8
)


# 11. Left panel: protein and localization ------------------

par(
  mar = c(6, 0.5, 3, 0),
  xaxs = "i",
  yaxs = "i"
)

plot(
  NA,
  xlim = c(0, 1),
  ylim = ylim,
  axes = FALSE,
  xlab = "",
  ylab = "",
  type = "n"
)

# Column positions
gene_x <- 0.06
location_x <- 0.34

# Font sizes
gene_cex <- 0.95
location_cex <- 0.80

# Ensure sufficient space for the longest protein name
gene_width <- max(
  strwidth(genes, cex = gene_cex)
)

location_x <- max(
  location_x,
  gene_x + gene_width + 0.035
)

# Ensure localization text fits horizontally
fit_label <- function(x, max_width, cex) {

  if (strwidth(x, cex = cex) <= max_width) {
    return(x)
  }

  # Remove secondary localization first
  x <- strsplit(x, " / ", fixed = TRUE)[[1]][1]

  if (strwidth(x, cex = cex) <= max_width) {
    return(x)
  }

  # Truncate only if necessary
  while (
    nchar(x) > 1 &&
    strwidth(
      paste0(x, "..."),
      cex = cex
    ) > max_width
  ) {
    x <- substr(x, 1, nchar(x) - 1)
  }

  paste0(x, "...")
}

location <- vapply(
  location,
  fit_label,
  max_width = 0.99 - location_x,
  cex = location_cex,
  FUN.VALUE = character(1)
)

# Protein names
text(
  x = gene_x,
  y = bp,
  labels = genes,
  adj = c(0, 0.5),
  font = 1,
  cex = gene_cex
)

# Localization annotations
text(
  x = location_x,
  y = bp,
  labels = location,
  adj = c(0, 0.5),
  col = "grey35",
  cex = location_cex
)

# Column headings
mtext(
  "Protein",
  side = 3,
  at = gene_x,
  adj = 0,
  line = 0.6,
  font = 2,
  cex = 0.85
)

mtext(
  "UniProt localization",
  side = 3,
  at = location_x,
  adj = 0,
  line = 0.6,
  font = 2,
  cex = 0.85
)


# 12. Right panel: fold-change barplot -----------------------

par(
  mar = c(6, 0, 3, 2),
  xaxs = "i",
  yaxs = "i"
)

bp <- as.numeric(
  barplot(
    values,
    names.arg = rep("", length(values)),
    horiz = TRUE,
    col = cols,
    density = ifelse(imputed, 25, NA),
    angle = 45,
    border = ifelse(imputed, cols, NA),
    xlim = xlim,
    ylim = ylim,
    axes = FALSE,
    axisnames = FALSE,
    width = 0.7,
    space = 0.65
  )
)

# Zero reference line
abline(
  v = 0,
  col = "grey35",
  lwd = 1
)

# Fold-change axis
axis(
  side = 1,
  at = seq(-axis_limit, axis_limit, 1),
  lwd = 1,
  lwd.ticks = 1,
  cex.axis = 0.85,
  tck = -0.015
)

# Axis label
mtext(
  expression(
    "Median-centered " * log[2] * "(nucleus/cytoplasm)"
  ),
  side = 1,
  line = 2.5,
  cex = 0.95
)


# 13. Numerical value labels --------------------------------

value_labels <- paste0(
  sprintf("%+.2f", values),
  ifelse(imputed, "*", "")
)

# Positive values: right of bars
pos <- values >= 0

text(
  x = values[pos],
  y = bp[pos],
  labels = value_labels[pos],
  pos = 4,
  offset = 0.35,
  cex = 0.8,
  xpd = NA
)

# Negative values: left of bars
neg <- values < 0

text(
  x = values[neg],
  y = bp[neg],
  labels = value_labels[neg],
  pos = 2,
  offset = 0.35,
  cex = 0.8,
  xpd = NA
)


# 14. Compartment labels -------------------------------------

mtext(
  "Cytoplasmic enrichment",
  side = 3,
  at = -plot_limit / 2,
  line = 0.6,
  col = col_cytoplasm,
  font = 2,
  cex = 0.85
)

mtext(
  "Nuclear enrichment",
  side = 3,
  at = plot_limit / 2,
  line = 0.6,
  col = col_nucleus,
  font = 2,
  cex = 0.85
)

# Imputation annotation
mtext(
  "* Missing intensity imputed",
  side = 1,
  line = 4.0,
  cex = 0.7,
  col = "grey40"
)


# 15. Save annotated data ------------------------------------

# Preserve the exact labels used in the figure
res$Localization_plot <- rev(location)

# Save complete UniProt annotations
write.table(
  uniprot,
  "./data/csv/LC-MSMS/UniProt_localization.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = TRUE
)

# Save annotated marker proteins
write.csv(
  res,
  "./data/csv/LC-MSMS/NORMA_annotated_markers.csv",
  row.names = FALSE
)

# Save complete normalized proteome
write.csv(
  dat,
  "./data/csv/LC-MSMS/NORMA_normalized_proteome.csv",
  row.names = FALSE
)

# Reset layout
layout(1)

quartz.save(file="./figures/pdf/uniprot_foldchange.pdf", type="pdf")
