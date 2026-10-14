#!/usr/bin/env Rscript
# Annotated cgMLST minimum spanning tree (GrapeTree MSTreeV2) + GrapeTree metadata file.
# Usage: Rscript bin/plot_mst.R <mst.nwk> <combined_summary.tsv> <clusters.tsv|NONE> <prefix>
# Output: <prefix>.pdf (one page per annotation) and grapetree_metadata.tsv
#         (load both cgmlst_mst.nwk and this file at https://achtman-lab.github.io/GrapeTree/).
# Failures are non-fatal: exits 0 and writes PLOT_MST_FAILED instead of crashing.

suppressMessages({
  library(ape)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggtree)
  library(RColorBrewer)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: plot_mst.R <mst.nwk> <summary.tsv> <clusters.tsv|NONE> <prefix>")
nwk <- args[1]; summary_file <- args[2]; cluster_file <- args[3]; prefix <- args[4]

top_n_lab <- function(x, n = 12) {           # top n categories keep their name; rest "Other"
  x <- ifelse(is.na(x) | x %in% c("", "-"), "Unknown", as.character(x))
  keep <- names(sort(table(x[x != "Unknown"]), decreasing = TRUE))[seq_len(n)]
  ifelse(x %in% c(keep, "Unknown"), x, "Other")
}

# ── Metadata: one row per sample, GrapeTree's first column must be "ID" ──────
s <- read.delim(summary_file, check.names = FALSE, quote = "", na.strings = c("", "NA"),
                colClasses = "character")
col <- function(name) if (name %in% names(s)) s[[name]] else NA_character_
meta <- data.frame(
  ID               = s$sample_id,
  ST               = col("mlst_st"),
  spa_type         = col("spa_type"),
  SCCmec           = col("sccmec_group"),
  agr_group        = col("agr_group"),
  PVL              = ifelse(str_detect(coalesce(col("virulence_summary"), ""), "PVL\\+"), "PVL+", "PVL-"),
  mecA             = ifelse(str_detect(coalesce(col("amrfinder_genes"), ""), "mecA"), "mecA+", "mecA-"),
  mecC             = ifelse(str_detect(coalesce(col("amrfinder_genes"), ""), "mecC"), "mecC+", "mecC-"),
  country          = col("geo_loc_country"),
  region           = col("geo_loc_region"),
  year             = col("year"),
  collection_date  = col("collection_date"),
  host             = col("host_group"),
  isolation_source = col("source_group"),
  stringsAsFactors = FALSE
)
if (file.exists(cluster_file)) {
  cl <- read.delim(cluster_file, check.names = FALSE, na.strings = c("", "NA"), colClasses = "character")
  meta <- left_join(meta, cl, by = c("ID" = "sample_id"))
}
write.table(meta, "grapetree_metadata.tsv", sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# ── MST plot ─────────────────────────────────────────────────────────────────
tryCatch({
  tree <- read.tree(nwk)
  alleles <- tree$edge.length
  # display length: log scale so non-CC lineages (1000+ alleles) don't squash the rest;
  # identical genomes (0 alleles) get a short stub so they fan out around their neighbour
  tree$edge.length <- log1p(alleles) + 0.4
  p0 <- ggtree(tree, layout = "equal_angle", colour = "grey55", linewidth = 0.3)
  d <- p0$data
  d$alleles <- c(alleles, NA)[match(d$node, c(tree$edge[, 2], NA))]
  par_xy <- d[match(d$parent, d$node), c("x", "y")]
  d$mx <- (d$x + par_xy$x) / 2
  d$my <- (d$y + par_xy$y) / 2
  tips <- left_join(filter(d, isTip), meta, by = c("label" = "ID"))

  pages <- list(
    "Country (top 12)"     = top_n_lab(tips$country),
    "ST (top 12)"          = top_n_lab(tips$ST),
    "SCCmec"               = top_n_lab(tips$SCCmec),
    "Host"                 = top_n_lab(tips$host),
    "Isolation source"     = top_n_lab(tips$isolation_source),
    "Collection year"      = ifelse(is.na(tips$year), "Unknown", tips$year),
    "PVL"                  = tips$PVL,
    "Outbreak cluster (top 12 multi-member)" = if ("outbreak_cluster" %in% names(tips)) {
      multi <- names(which(table(tips$outbreak_cluster) > 1))
      sub("^Unknown$", "Not clustered", top_n_lab(ifelse(tips$outbreak_cluster %in% multi, tips$outbreak_cluster, NA)))
    } else rep("Unknown", nrow(tips))
  )

  edge_lab <- filter(d, !is.na(alleles), alleles > 0)
  pdf(paste0(prefix, ".pdf"), width = 14, height = 14)
  for (title in names(pages)) {
    tips$var <- pages[[title]]
    lv <- sort(setdiff(unique(tips$var), c("Other", "Unknown")))
    pal <- if (title == "Collection year") colorRampPalette(c("#fff5b1", "#fd8d3c", "#7f0000"))(max(1, length(lv)))
           else colorRampPalette(c(brewer.pal(9, "Set1")[-c(6, 9)], brewer.pal(8, "Dark2")))(max(1, length(lv)))   # no yellow/grey
    lv <- setdiff(lv, "Not clustered")
    cols <- c(setNames(pal[seq_along(lv)], lv), Other = "#969696", Unknown = "#f0f0f0", "Not clustered" = "#ffffff")
    print(p0 +
      geom_text(data = edge_lab, aes(x = mx, y = my, label = alleles), size = 1.6, colour = "grey30") +
      geom_point(data = tips, aes(x = x, y = y, fill = var), shape = 21, size = 2.4, stroke = 0.2) +
      scale_fill_manual(values = cols, breaks = c(lv, "Other", "Unknown", "Not clustered"), name = NULL) +
      labs(title = sprintf("cgMLST minimum spanning tree (GrapeTree MSTreeV2, %d genomes): %s", nrow(tips), title),
           caption = "Numbers on edges = allele differences. Edge lengths drawn on a log scale; identical genomes sit on short stubs.") +
      theme(legend.position = "right", plot.title = element_text(size = 11, face = "bold")))
  }
  invisible(dev.off())
  message("Saved: ", prefix, ".pdf (", length(pages), " pages)")
}, error = function(e) {
  message("MST plot failed: ", conditionMessage(e))
  writeLines(conditionMessage(e), "PLOT_MST_FAILED")
})
