#!/usr/bin/env Rscript
# summarize_coverage_quantized.R
#
# For each quantized BED file in results/mosdepth/ (auto-discovered):
#   1. Load <strain>.quantized.bed.gz (BED4: chrom, start, end, label)
#      Labels: NO_COVERAGE, LOW_COVERAGE, CALLABLE, HIGH_COVERAGE, VERY_HIGH_COVERAGE
#      Covers plain Illumina, _ont, and _pb variants automatically.
#   2. Produce per-strain multi-page PDFs (50 scaffolds per page):
#        heatmap - geom_rect tiles coloured by coverage class, all scaffolds stacked
#        bar     - stacked bar showing proportion of each class per scaffold
#
# Output: results/mosdepth_quantized/<strain>/{heatmap,bar}.pdf

library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(forcats)
library(Biostrings)
library(parallel)

# ---- paths ----
mosdepth_dir  <- file.path("results", "mosdepth")
genome_dir    <- file.path("genome")
out_base      <- file.path("results", "mosdepth_quantized")
dir.create(out_base, showWarnings = FALSE, recursive = TRUE)

SCAFFOLDS_PER_PAGE <- 50
MIN_SCAFFOLD_LENGTH <- 50000
GC_WINDOW <- 10000   # bp window for GC% smoothing

# ---- discover all quantized BED files (plain, _ont, _pb variants) ----
bed_files <- list.files(mosdepth_dir, pattern = "\\.quantized\\.bed\\.gz$",
                        full.names = FALSE)
if (length(bed_files) == 0) stop("No quantized BED files found in ", mosdepth_dir)
strains <- sub("\\.quantized\\.bed\\.gz$", "", bed_files)
message("Found ", length(strains), " quantized BED files: ",
        paste(strains, collapse = ", "))

# ---- quantized coverage levels: cool -> hot ----
cov_levels <- c("NO_COVERAGE", "LOW_COVERAGE", "CALLABLE",
                "HIGH_COVERAGE", "VERY_HIGH_COVERAGE")

# cool-to-hot palette matching the five classes
cov_colours <- c(
  NO_COVERAGE        = "#313695",   # deep blue  (cool)
  LOW_COVERAGE       = "#74add1",   # light blue
  CALLABLE           = "#ffffbf",   # pale yellow (neutral)
  HIGH_COVERAGE      = "#f46d43",   # orange
  VERY_HIGH_COVERAGE = "#a50026"    # deep red   (hot)
)

# ---- shared theme ----
base_theme <- theme_bw(base_size = 10) +
  theme(
    strip.text       = element_text(face = "bold", size = 7),
    strip.background = element_rect(fill = "#e8eaf6"),
    panel.spacing    = unit(0.3, "lines"),
    legend.position  = "right",
    axis.title       = element_text(size = 9)
  )

# ---- helper: split a vector into chunks of n ----
chunk <- function(x, n) split(x, ceiling(seq_along(x) / n))

# ---- helper: compute GC% in non-overlapping windows (vectorized via Views) ----
compute_gc_windows <- function(genome, scaffolds, window_size = GC_WINDOW) {
  scaffolds_present <- intersect(scaffolds, names(genome))
  if (length(scaffolds_present) == 0) return(NULL)
  result <- lapply(scaffolds_present, function(chr) {
    sq     <- genome[[chr]]
    len    <- length(sq)
    starts <- seq(1, len, by = window_size)
    ends   <- pmin(starts + window_size - 1, len)
    views  <- Views(sq, start = starts, end = ends)
    gc_pct <- letterFrequency(views, letters = "GC", as.prob = TRUE)[, 1] * 100
    data.frame(chrom = chr,
               pos   = (starts + ends) / 2 - 1,  # midpoint, 0-based
               gc    = gc_pct)
  })
  bind_rows(result)
}

# ---- helper: resolve genome FASTA for a strain (strip _ont / _pb suffixes) ----
find_genome_fasta <- function(strain) {
  base <- sub("_(ont|pb)$", "", strain)
  path <- file.path(genome_dir, paste0(base, ".masked.fasta"))
  if (file.exists(path)) return(path)
  NULL
}

# ===========================================================================
# Per-strain processing function (parallelised below)
# ===========================================================================
process_strain <- function(strain) {

  message("\n=== ", strain, " ===")

  # ---- load quantized BED ----
  bed_file <- file.path(mosdepth_dir, paste0(strain, ".quantized.bed.gz"))
  if (!file.exists(bed_file)) {
    warning("Quantized BED file not found, skipping: ", bed_file)
    return(invisible(NULL))
  }

  bed <- read_tsv(bed_file,
                  col_names = c("chrom", "start", "end", "label"),
                  col_types = cols(
                    chrom = col_character(),
                    start = col_double(),
                    end   = col_double(),
                    label = col_character()
                  ))

  # ---- set factor levels; scaffold order = order of appearance in BED ----
  scaffold_order <- unique(bed$chrom)

  # ---- filter scaffolds shorter than MIN_SCAFFOLD_LENGTH ----
  scaffold_lengths <- bed |>
    group_by(chrom) |>
    summarise(length = max(end) - min(start), .groups = "drop")
  keep_scaffolds <- scaffold_lengths |>
    filter(length >= MIN_SCAFFOLD_LENGTH) |>
    pull(chrom)
  n_dropped      <- length(scaffold_order) - length(keep_scaffolds)
  scaffold_order <- scaffold_order[scaffold_order %in% keep_scaffolds]
  if (n_dropped > 0)
    message("  dropped ", n_dropped, " scaffolds shorter than ", MIN_SCAFFOLD_LENGTH, " bp")

  if (length(scaffold_order) == 0) {
    warning("No scaffolds remain after length filter for ", strain)
    return(invisible(NULL))
  }

  dat <- bed |>
    filter(chrom %in% scaffold_order) |>
    mutate(
      chrom = factor(chrom, levels = scaffold_order),
      label = factor(label, levels = cov_levels),
      width = end - start
    )

  nscaff      <- length(scaffold_order)
  out_dir     <- file.path(out_base, strain)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  strain_label <- str_replace_all(strain, "\\.", " ")
  pages        <- chunk(scaffold_order, SCAFFOLDS_PER_PAGE)
  npages       <- length(pages)

  # ---- load genome and compute GC% windows ----
  gc_dat <- NULL
  genome_fasta <- find_genome_fasta(strain)
  if (!is.null(genome_fasta)) {
    message("  computing GC% from ", basename(genome_fasta))
    genome <- readDNAStringSet(genome_fasta)
    names(genome) <- sub(" .*", "", names(genome))
    gc_dat <- compute_gc_windows(genome, scaffold_order)
    if (!is.null(gc_dat))
      gc_dat$chrom <- factor(gc_dat$chrom, levels = scaffold_order)
  } else {
    message("  no genome FASTA found for ", strain, " — GC panel will be skipped")
  }

  message("  ", nscaff, " scaffolds -> ", npages, " pages (", SCAFFOLDS_PER_PAGE, " per page)")

  # ---- pre-split page data once to avoid repeated full-table scans ----
  page_dat_list <- lapply(pages, function(scaffs)
    dat |> filter(chrom %in% scaffs) |> mutate(chrom = factor(chrom, levels = scaffs)))

  page_gc_list <- if (!is.null(gc_dat))
    lapply(pages, function(scaffs)
      gc_dat |> filter(chrom %in% scaffs) |> mutate(chrom = factor(chrom, levels = scaffs)))
  else
    vector("list", npages)

  # =========================================================================
  # 1. HEATMAP  (geom_rect tiles, one row per scaffold, colour = class)
  #    Multi-page PDF: 50 scaffolds per page
  # =========================================================================
  heat_path <- file.path(out_dir, "heatmap.pdf")
  cairo_pdf(heat_path, width = 14, height = 20, onefile = TRUE)

  for (i in seq_along(pages)) {
    page_dat <- page_dat_list[[i]]
    page_gc  <- page_gc_list[[i]]

    # Base heatmap — GC overlay added below when available
    p_heat <- ggplot(page_dat, aes(xmin = start, xmax = end,
                                   ymin = 0,    ymax = 1,
                                   fill = label)) +
      geom_rect() +
      scale_fill_manual(values = cov_colours,
                        drop   = FALSE,
                        name   = "Coverage\nclass") +
      facet_wrap(~chrom, ncol = 1, scales = "free_x",
                 strip.position = "left") +
      labs(
        title    = paste0(strain_label,
                          " - scaffold coverage (quantized heatmap)",
                          "  [page ", i, " / ", npages, "]"),
        subtitle = "Colour = mosdepth quantized coverage class",
        x        = "Genomic position (bp)",
        y        = NULL
      ) +
      base_theme +
      theme(
        strip.text.y.left  = element_text(angle = 0, hjust = 1),
        axis.text.y.left   = element_blank(),
        axis.ticks.y.left  = element_blank(),
        panel.grid         = element_blank()
      )

    if (!is.null(page_gc)) {
      gc_min    <- min(page_gc$gc, na.rm = TRUE)
      gc_max    <- max(page_gc$gc, na.rm = TRUE)
      gc_median <- median(page_gc$gc, na.rm = TRUE)
      gc_range  <- gc_max - gc_min
      if (gc_range < 0.01) gc_range <- 1  # guard against flat GC (e.g. all-N scaffolds)
      gc_breaks <- pretty(c(gc_min, gc_max), n = 4)
      gc_breaks <- gc_breaks[gc_breaks >= gc_min & gc_breaks <= gc_max]

      page_gc_scaled <- page_gc |>
        mutate(gc_scaled = (gc - gc_min) / gc_range)

      p_heat <- p_heat +
        geom_hline(yintercept = (gc_median - gc_min) / gc_range,
                   colour = "grey60", linewidth = 0.3, linetype = "dashed") +
        geom_line(data        = page_gc_scaled,
                  mapping     = aes(x = pos, y = gc_scaled),
                  inherit.aes = FALSE,
                  colour      = "#2ca25f", linewidth = 0.4) +
        scale_y_continuous(
          limits   = c(0, 1),
          sec.axis = sec_axis(~ . * gc_range + gc_min, name = "GC%",
                              breaks = (gc_breaks - gc_min) / gc_range,
                              labels = paste0(round(gc_breaks, 1), "%"))
        ) +
        labs(subtitle = "Colour = mosdepth quantized coverage class; green line = GC%")
    }

    print(p_heat)
  }

  dev.off()
  message("  saved -> ", heat_path)

  # =========================================================================
  # 2. STACKED BAR  (proportion of each class per scaffold)
  #    Multi-page PDF: 50 scaffolds per page
  # =========================================================================
  bar_dat <- dat |>
    group_by(chrom, label) |>
    summarise(total_bp = sum(width), .groups = "drop") |>
    group_by(chrom) |>
    mutate(pct = total_bp / sum(total_bp) * 100) |>
    ungroup()

  bar_path <- file.path(out_dir, "bar.pdf")
  cairo_pdf(bar_path, width = 12, height = 18, onefile = TRUE)

  for (i in seq_along(pages)) {
    scaffs       <- pages[[i]]
    page_bar_dat <- bar_dat |> filter(chrom %in% scaffs) |>
                      mutate(chrom = factor(chrom, levels = rev(scaffs)))

    p_bar <- ggplot(page_bar_dat, aes(x = pct, y = chrom, fill = label)) +
      geom_col(width = 0.8) +
      scale_fill_manual(values = cov_colours,
                        drop   = FALSE,
                        name   = "Coverage\nclass") +
      scale_x_continuous(expand = c(0, 0),
                         labels = function(x) paste0(x, "%")) +
      labs(
        title    = paste0(strain_label,
                          " - scaffold coverage class proportions",
                          "  [page ", i, " / ", npages, "]"),
        subtitle = "Percentage of scaffold bases in each quantized class",
        x        = "Percentage of scaffold (bp)",
        y        = "Scaffold"
      ) +
      base_theme

    print(p_bar)
  }

  dev.off()
  message("  saved -> ", bar_path)

  message("  done: ", nscaff, " scaffolds plotted across ", npages, " pages")
  invisible(NULL)
}

# ===========================================================================
# Run: parallel over strains (falls back to sequential if only one core)
# ===========================================================================
ncores <- min(length(strains), max(1L, detectCores() - 1L))
message("Processing ", length(strains), " strain(s) on ", ncores, " core(s)")

if (ncores > 1) {
  mclapply(strains, process_strain, mc.cores = ncores)
} else {
  lapply(strains, process_strain)
}

message("\nAll strains processed. Output in: ", out_base)
