#!/usr/bin/env Rscript
# Annotated core-SNP phylogenetic tree from staphit2 results.
# Usage: Rscript bin/plot_tree.R <results_subdir> [run_name]
# Output: <run_name>_tree.pdf  and  <run_name>_tree.svg  in the current directory.
# Failures are non-fatal: exits 0 and writes PLOT_TREE_FAILED instead of crashing.

# ── Dependencies ─────────────────────────────────────────────────────────────
local_lib <- path.expand("~/.R/library")
dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_lib, .libPaths()))

install_if_missing <- function(pkgs, bioc = FALSE) {
  missing <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
  if (!length(missing)) return(invisible())
  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repos = "https://cloud.r-project.org",
                       lib = local_lib, quiet = TRUE)
    BiocManager::install(missing, ask = FALSE, quiet = TRUE, lib = local_lib)
  } else {
    install.packages(missing, repos = "https://cloud.r-project.org",
                     lib = local_lib, quiet = TRUE)
  }
}

# Wrap everything so R package issues never crash a Nextflow pipeline
tryCatch({

install_if_missing(c("dplyr", "readr", "tidyr", "stringr", "ggplot2",
                     "ggnewscale", "RColorBrewer", "svglite"))
install_if_missing(c("treeio", "ggtree", "ggtreeExtra"), bioc = TRUE)

suppressMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(treeio)
  library(ggtree)
  library(ggtreeExtra)
  library(ggnewscale)
  library(RColorBrewer)
  library(svglite)
})

# ── Arguments ─────────────────────────────────────────────────────────────────
args     <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Usage: Rscript bin/plot_tree.R <results_subdir> [run_name]")
outdir   <- args[1]
run_name <- if (length(args) >= 2) args[2] else basename(outdir)

tree_file    <- file.path(outdir, "iqtree",   "core.treefile")
summary_file <- file.path(outdir, "summary",  "combined_summary.tsv")
cluster_file <- file.path(outdir, "clusters", "clusters.tsv")

if (!file.exists(tree_file))    stop("Tree not found: ", tree_file)
if (!file.exists(summary_file)) stop("Summary not found: ", summary_file)

# ── Load & parse metadata ──────────────────────────────────────────────────────
tree <- read.tree(tree_file)

meta <- read_tsv(summary_file, show_col_types = FALSE) |>
  mutate(
    pvl = case_when(
      str_detect(virulence_summary, "PVL\\+") ~ "PVL+",
      str_detect(virulence_summary, "PVL-")   ~ "PVL-",
      TRUE ~ NA_character_
    ),
    tsst = case_when(
      str_detect(virulence_summary, "TSST\\+") ~ "TSST+",
      str_detect(virulence_summary, "TSST-")   ~ "TSST-",
      TRUE ~ NA_character_
    ),
    meca = if_else(
      str_detect(coalesce(amrfinder_genes, ""), "mecA"),
      "mecA+", "mecA-"
    ),
    sccmec = case_when(
      is.na(sccmec_type) | sccmec_type %in% c("Negative", "ND") ~ "MSSA",
      TRUE ~ sccmec_type
    ),
    spa_group = {
      top12 <- names(sort(table(spa_type), decreasing = TRUE))[1:12]
      if_else(spa_type %in% top12, spa_type, "Other")
    }
  )

if (file.exists(cluster_file)) {
  meta <- left_join(meta,
    read_tsv(cluster_file, show_col_types = FALSE) |>
      select(sample_id, outbreak_cluster),
    by = "sample_id"
  )
} else {
  meta$outbreak_cluster <- NA_character_
  message("No clusters.tsv found — outbreak strip will be empty")
}

multi_ob <- meta |>
  filter(!is.na(outbreak_cluster)) |>
  count(outbreak_cluster) |>
  filter(n > 1) |>
  pull(outbreak_cluster)

meta <- meta |>
  mutate(ob_strip = if_else(outbreak_cluster %in% multi_ob, outbreak_cluster, NA_character_))

# ── Colour palettes ───────────────────────────────────────────────────────────
spa_levels  <- setdiff(unique(na.omit(meta$spa_group)), "Other")
spa_palette <- colorRampPalette(brewer.pal(8, "Set2"))(length(spa_levels))
spa_colors  <- c(setNames(spa_palette, spa_levels), Other = "grey75")

sccmec_colors <- c(
  MSSA                  = "#d9d9d9",
  "Type IV"             = "#4daf4a",
  "Type V"              = "#377eb8",
  "Type VII"            = "#ff7f00",
  "Type XI"             = "#984ea3",
  "Composite (Type IV)" = "#a8d978",
  Unknown               = "#eeeeee"
)

agr_colors  <- c(gp1 = "#1b9e77", gp2 = "#d95f02", gp3 = "#7570b3", gp4 = "#e7298a")
pvl_colors  <- c("PVL+"  = "#d73027", "PVL-"  = "#f5f5f5")
tsst_colors <- c("TSST+" = "#4575b4", "TSST-" = "#f5f5f5")
meca_colors <- c("mecA+" = "#b2182b", "mecA-" = "#f5f5f5")

# Only keep top 12 outbreak clusters by size to avoid legend explosion
ob_top12 <- meta |>
  filter(!is.na(outbreak_cluster)) |>
  count(outbreak_cluster, sort = TRUE) |>
  filter(n > 1) |>
  slice_head(n = 12) |>
  pull(outbreak_cluster)

meta <- meta |>
  mutate(ob_strip = if_else(outbreak_cluster %in% ob_top12, outbreak_cluster, NA_character_))

n_ob <- length(ob_top12)
ob_colors <- if (n_ob == 0) {
  character(0)
} else {
  pal <- brewer.pal(max(3, min(n_ob, 12)), "Paired")
  setNames(pal[seq_len(n_ob)], sort(ob_top12))
}

# ── Helper: add one annotation ring (circular layout) ────────────────────────
add_ring <- function(p, col, fill_scale, offset = 0.05, pwidth = 0.1) {
  df <- meta |> select(sample_id, strip_val = all_of(col))
  p + new_scale_fill() +
    geom_fruit(
      data    = df,
      geom    = geom_tile,
      mapping = aes(y = sample_id, fill = strip_val),
      offset  = offset,
      pwidth  = pwidth,
      color   = NA
    ) +
    fill_scale
}

# ── Build circular tree ───────────────────────────────────────────────────────
p <- ggtree(tree, layout = "circular", linewidth = 0.15, color = "grey40") %<+% meta +
  geom_tippoint(aes(color = spa_group), size = 0.6, alpha = 0.9) +
  scale_color_manual(values = spa_colors, name = "spa type", na.value = "grey80") +
  theme_tree() +
  theme(
    legend.position  = "right",
    legend.key.size  = unit(0.3, "cm"),
    legend.text      = element_text(size = 6.5),
    legend.title     = element_text(size = 7.5, face = "bold"),
    legend.spacing.y = unit(0.1, "cm"),
    plot.title       = element_text(size = 10, face = "bold", hjust = 0.5),
    plot.caption     = element_text(size = 6, color = "grey55", hjust = 0.5),
    plot.margin      = margin(10, 10, 10, 10)
  ) +
  labs(
    title   = sprintf("%s - Core SNP phylogeny (%d isolates)", run_name, Ntip(tree)),
    caption = "Rings (inner->outer): SCCmec | agr group | PVL | TSST | mecA | Top-12 outbreak clusters"
  )

p <- add_ring(p, "sccmec",
  scale_fill_manual(values = sccmec_colors, name = "SCCmec",    na.value = "#eeeeee"))
p <- add_ring(p, "agr_group",
  scale_fill_manual(values = agr_colors,    name = "agr group", na.value = "grey90"))
p <- add_ring(p, "pvl",
  scale_fill_manual(values = pvl_colors,    name = "PVL",       na.value = "grey90"))
p <- add_ring(p, "tsst",
  scale_fill_manual(values = tsst_colors,   name = "TSST",      na.value = "grey90"))
p <- add_ring(p, "meca",
  scale_fill_manual(values = meca_colors,   name = "mecA",      na.value = "grey90"))

if (n_ob > 0) {
  p <- add_ring(p, "ob_strip",
    scale_fill_manual(values = ob_colors, name = "Outbreak\ncluster", na.value = "grey97"))
}

# ── Save ──────────────────────────────────────────────────────────────────────
width  <- 14
height <- 14

out_pdf <- paste0(run_name, "_tree.pdf")
out_svg <- paste0(run_name, "_tree.svg")

ggsave(out_pdf, p, width = width, height = height, limitsize = FALSE)
ggsave(out_svg, p, width = width, height = height, limitsize = FALSE,
       device = svglite, fix_text_size = FALSE)

message(sprintf("Saved: %s  (%.0f x %.0f in)", out_pdf, width, height))
message(sprintf("Saved: %s  (%.0f x %.0f in)", out_svg, width, height))

}, error = function(e) {
  msg <- conditionMessage(e)
  message("WARNING: plot_tree.R failed — ", msg)
  message("Tree visualisation skipped. Pipeline output is unaffected.")
  writeLines(msg, "PLOT_TREE_FAILED")
  quit(status = 0)
})
