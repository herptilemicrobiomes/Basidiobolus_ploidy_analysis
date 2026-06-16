#!/usr/bin/env Rscript
# Chromosome/contig-organized per-gene ploidy plots from nQuire per-gene output.
#
# For each strain it joins results/nquire2/<strain>.pergene.nquire.tsv.gz to
# bed/<strain>.bed (gene -> contig,start,end), picks the winning ploidy model as
# the SMALLEST delta-log-likelihood among d_dip / d_tri / d_tet (0 and nan are
# treated as missing/unestimable), and writes a 2-page PDF per strain:
#   page 1 = top-N contigs by size, faceted, genes ordered by position
#   page 2 = genome-wide Manhattan track, contigs ordered by size,
#            alternating red/black so genes on the same contig group together
#
# Usage:
#   Rscript pipeline/plot/plot_pergene_ploidy.R [--mode ploidy|score]
#                                               [--top 50] [--strain NAME]...
#   --mode ploidy  (default) y-axis = winning ploidy call (2n/3n/4n)
#   --mode score             y-axis = winning delta-log-likelihood score
#
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

## ---------------------------------------------------------------- args ----
args   <- commandArgs(trailingOnly = TRUE)
mode   <- "ploidy"
top_n  <- 50L
strains_req <- NULL
i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--mode")        { mode <- args[i + 1L];                 i <- i + 2L }
  else if (a == "--top")    { top_n <- as.integer(args[i + 1L]);    i <- i + 2L }
  else if (a == "--strain") { strains_req <- c(strains_req, args[i + 1L]); i <- i + 2L }
  else                      { i <- i + 1L }
}
stopifnot(mode %in% c("ploidy", "score"))

bed_dir <- "bed"
nq_dir  <- "results/nquire2"
fai_dir <- "genome"
out_dir <- "results/nquire2/plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## --------------------------------------------------------- definitions ----
models <- c(d_dip = "2n (diploid)", d_tri = "3n (triploid)", d_tet = "4n (tetraploid)")
plevel <- c(d_dip = 2, d_tri = 3, d_tet = 4)
ploidy_cols <- c("2n (diploid)"    = "#1b9e77",
                 "3n (triploid)"   = "#d95f02",
                 "4n (tetraploid)" = "#7570b3")
ploidy_levels <- unname(models)

## ----------------------------------------------------------- per strain ----
files <- list.files(nq_dir, pattern = "\\.pergene\\.nquire\\.tsv\\.gz$", full.names = TRUE)
all_strains <- sub("\\.pergene\\.nquire\\.tsv\\.gz$", "", basename(files))
if (!is.null(strains_req)) all_strains <- intersect(all_strains, strains_req)

summary_rows  <- list()
manhattan_all <- list()

for (strain in all_strains) {
  nqf  <- file.path(nq_dir,  paste0(strain, ".pergene.nquire.tsv.gz"))
  bedf <- file.path(bed_dir, paste0(strain, ".bed"))
  faif <- file.path(fai_dir, paste0(strain, ".masked.fasta.fai"))
  if (!file.exists(bedf)) { message("skip ", strain, ": no bed file"); next }

  dt <- tryCatch(fread(nqf, na.strings = c("nan", "-nan", "NA", "")),
                 error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0L) { message("skip ", strain, ": empty per-gene table"); next }
  bed <- fread(bedf, header = FALSE, col.names = c("contig", "start", "end", "gene"))

  ## ---- winner = smallest delta-LL; 0 and non-finite are missing ----------
  m <- as.matrix(dt[, c("d_dip", "d_tri", "d_tet")])
  m[!is.finite(m)] <- NA
  m[m == 0] <- NA
  all_na  <- rowSums(!is.na(m)) == 0L
  wins    <- max.col(-replace(m, is.na(m), Inf), ties.method = "first")  # min delta-LL
  win_mod <- colnames(m)[wins]
  win_scr <- m[cbind(seq_len(nrow(m)), wins)]
  win_mod[all_na] <- NA
  win_scr[all_na] <- NA

  res <- data.table(
    gene       = dt$gene,
    ploidy_lab = factor(models[win_mod], levels = ploidy_levels),
    score      = win_scr
  )
  res <- merge(bed, res, by = "gene")
  if (nrow(res) == 0L) { message("skip ", strain, ": no gene-name overlap bed<->nquire"); next }

  ## ---- contig lengths (fai if present, else max gene end) ----------------
  if (file.exists(faif)) {
    ctab <- fread(faif, header = FALSE, select = 1:2, col.names = c("contig", "clen"))
  } else {
    ctab <- bed[, .(clen = max(end)), by = contig]
  }
  setorder(ctab, -clen)
  ctab[, rank := .I]

  res <- merge(res, ctab, by = "contig")
  res[, called := !is.na(ploidy_lab)]

  ## ---- console summary ---------------------------------------------------
  tab <- res[, .N, by = ploidy_lab][order(ploidy_lab)]
  message(sprintf("== %s : %d genes, %d called, %d missing ==",
                  strain, nrow(res), sum(res$called), sum(!res$called)))
  print(tab)
  summary_rows[[strain]] <- data.table(
    strain  = strain,
    total   = nrow(res),
    called  = sum(res$called),
    missing = sum(!res$called),
    n2 = sum(res$ploidy_lab == "2n (diploid)",    na.rm = TRUE),
    n3 = sum(res$ploidy_lab == "3n (triploid)",   na.rm = TRUE),
    n4 = sum(res$ploidy_lab == "4n (tetraploid)", na.rm = TRUE)
  )

  is_score <- (mode == "score")

  ## ================= PAGE 1 : top-N contigs faceted ======================
  top_contigs <- ctab[rank <= top_n, contig]
  d1 <- res[contig %in% top_contigs & called == TRUE]
  ## facet label: contig (length kb, n genes), ordered by size
  lab_dt <- res[contig %in% top_contigs, .(ng = .N), by = contig]
  lab_dt <- merge(lab_dt, ctab, by = "contig")
  setorder(lab_dt, rank)
  lab_dt[, lab := sprintf("%s  (%.0f kb, %d genes)", contig, clen / 1e3, ng)]
  d1[, contig_f := factor(contig, levels = lab_dt$contig, labels = lab_dt$lab)]
  d1[, ploidy_lvl := plevel[c("2n (diploid)" = "d_dip", "3n (triploid)" = "d_tri",
                              "4n (tetraploid)" = "d_tet")[as.character(ploidy_lab)]]]
  d1[, posk := start / 1e3]

  if (is_score) {
    p1 <- ggplot(d1, aes(posk, score, colour = ploidy_lab)) +
      geom_point(size = 0.7, alpha = 0.8) +
      labs(y = "winning delta-log-likelihood (smaller = better fit)")
  } else {
    p1 <- ggplot(d1, aes(posk, ploidy_lvl, colour = ploidy_lab)) +
      geom_point(size = 0.8, alpha = 0.8,
                 position = position_jitter(height = 0.12, width = 0)) +
      scale_y_continuous(breaks = c(2, 3, 4),
                         labels = c("2n", "3n", "4n"),
                         limits = c(1.6, 4.4)) +
      labs(y = "winning ploidy call")
  }
  p1 <- p1 +
    facet_wrap(~contig_f, scales = "free_x") +
    scale_colour_manual(values = ploidy_cols, drop = FALSE, name = "ploidy call") +
    labs(x = "position on contig (kb)",
         title = sprintf("%s  -  per-gene ploidy, top %d contigs by size", strain, top_n),
         subtitle = "winner = smallest delta-log-likelihood among d_dip/d_tri/d_tet (0/nan = missing, omitted)") +
    theme_bw(base_size = 8) +
    theme(legend.position = "top",
          panel.grid.minor = element_blank(),
          strip.text = element_text(size = 6))

  ## ================= PAGE 2 : genome-wide Manhattan ======================
  d2 <- res[called == TRUE]
  ## cumulative x position, contigs ordered by size
  off <- ctab[, .(contig, clen, rank)]
  setorder(off, rank)
  off[, cumstart := cumsum(shift(clen, fill = 0))]
  d2 <- merge(d2, off[, .(contig, cumstart, rank)], by = "contig",
              suffixes = c("", ".off"))
  d2[, gx := (cumstart + start) / 1e6]                 # Mb
  d2[, band := factor(rank %% 2L)]                      # alternating red/black
  d2[, ploidy_lvl := plevel[c("2n (diploid)" = "d_dip", "3n (triploid)" = "d_tri",
                              "4n (tetraploid)" = "d_tet")[as.character(ploidy_lab)]]]
  setorder(d2, rank, start)
  manhattan_all[[strain]] <- d2[, .(strain = strain, contig, gx, band, ploidy_lvl, score)]

  if (is_score) {
    p2 <- ggplot(d2, aes(gx, score, colour = band)) +
      geom_point(size = 0.45, alpha = 0.7) +
      labs(y = "winning delta-log-likelihood")
  } else {
    p2 <- ggplot(d2, aes(gx, ploidy_lvl, colour = band)) +
      geom_point(size = 0.45, alpha = 0.7,
                 position = position_jitter(height = 0.12, width = 0)) +
      scale_y_continuous(breaks = c(2, 3, 4), labels = c("2n", "3n", "4n"),
                         limits = c(1.6, 4.4)) +
      labs(y = "winning ploidy call")
  }
  p2 <- p2 +
    scale_colour_manual(values = c("0" = "black", "1" = "#cc2a2a"), guide = "none") +
    labs(x = "genome position (Mb)  -  contigs ordered by size, left = largest",
         title = sprintf("%s  -  genome-wide per-gene ploidy (Manhattan)", strain),
         subtitle = "alternating black/red marks consecutive contigs; winner = smallest delta-log-likelihood") +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())

  ## ----------------------------------------------------------- write ------
  outf <- file.path(out_dir, sprintf("%s.pergene_ploidy.%s.pdf", strain, mode))
  pdf(outf, width = 16, height = 10, onefile = TRUE)
  print(p1)
  print(p2)
  invisible(dev.off())
  message("wrote ", outf)
}

if (length(summary_rows)) {
  sm <- rbindlist(summary_rows)
  message("\n==================== ploidy call summary (", mode, ") ====================")
  print(sm)
  fwrite(sm, file.path(out_dir, "ploidy_call_summary.tsv"), sep = "\t")
}

## ================== COMBINED multi-strain overview PDF ====================
if (length(manhattan_all)) {
  mha <- rbindlist(manhattan_all)
  ## strains ordered by overall fraction tetraploid (most polyploid on top)
  ord <- sm[, frac4 := n4 / pmax(called, 1)][order(-frac4), strain]
  mha[, strain := factor(strain, levels = ord)]

  ## ---- page 1: stacked proportion of ploidy calls per strain ------------
  pm <- melt(sm, id.vars = "strain", measure.vars = c("n2", "n3", "n4"),
             variable.name = "k", value.name = "n")
  pm[, ploidy_lab := factor(c(n2 = "2n (diploid)", n3 = "3n (triploid)",
                              n4 = "4n (tetraploid)")[as.character(k)],
                            levels = ploidy_levels)]
  pm[, strain := factor(strain, levels = ord)]
  pbar <- ggplot(pm, aes(strain, n, fill = ploidy_lab)) +
    geom_col(position = "fill") +
    scale_fill_manual(values = ploidy_cols, name = "ploidy call") +
    scale_y_continuous(labels = percent_format()) +
    labs(x = NULL, y = "fraction of estimable genes",
         title = "Per-gene ploidy composition across strains",
         subtitle = "winner = smallest delta-log-likelihood (d_dip/d_tri/d_tet); missing genes excluded") +
    coord_flip() +
    theme_bw(base_size = 12) +
    theme(legend.position = "top")

  ## ---- page 2: genome-wide ploidy, faceted by strain --------------------
  if (mode == "score") {
    pcmp <- ggplot(mha, aes(gx, score, colour = band)) +
      geom_point(size = 0.35, alpha = 0.6) +
      labs(y = "winning delta-log-likelihood")
  } else {
    pcmp <- ggplot(mha, aes(gx, ploidy_lvl, colour = band)) +
      geom_point(size = 0.35, alpha = 0.6,
                 position = position_jitter(height = 0.12, width = 0)) +
      scale_y_continuous(breaks = c(2, 3, 4), labels = c("2n", "3n", "4n"),
                         limits = c(1.6, 4.4)) +
      labs(y = "winning ploidy call")
  }
  pcmp <- pcmp +
    facet_grid(strain ~ ., scales = "free_x", switch = "y") +
    scale_colour_manual(values = c("0" = "black", "1" = "#cc2a2a"), guide = "none") +
    labs(x = "genome position (Mb)  -  contigs ordered by size within each strain",
         title = "Genome-wide per-gene ploidy across strains (Manhattan)",
         subtitle = "alternating black/red marks consecutive contigs") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.text.y.left = element_text(angle = 0))

  outc <- file.path(out_dir, sprintf("ALL_strains.pergene_ploidy.%s.pdf", mode))
  pdf(outc, width = 16, height = 11, onefile = TRUE)
  print(pbar)
  print(pcmp)
  invisible(dev.off())
  message("wrote ", outc)
}
