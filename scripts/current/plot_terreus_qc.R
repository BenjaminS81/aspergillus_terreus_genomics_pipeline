#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: plot_terreus_qc.R SUMMARY.tsv OUTPUT_DIR")

input_file <- args[[1]]
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(scales)
})

d <- read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE)
numeric_cols <- setdiff(names(d), c("label", "species_group", "sample", "source_fasta"))
for (col in numeric_cols) d[[col]] <- as.numeric(d[[col]])

species_levels <- sort(unique(d$species_group))
palette_values <- c(
  alabamensis = "#0072B2", aureoterreus = "#E69F00", floccosus = "#009E73",
  hortae = "#CC79A7", pseudoterreus = "#56B4E9", terreus = "#D55E00"
)
missing_species <- setdiff(species_levels, names(palette_values))
if (length(missing_species)) {
  palette_values <- c(palette_values, setNames(hue_pal()(length(missing_species)), missing_species))
}

d$display_label <- paste(d$species_group, d$sample, sep = " / ")
order_idx <- order(d$species_group, d$busco_complete_percent, d$sample)
d$display_label <- factor(d$display_label, levels = d$display_label[order_idx])
d$species_group <- factor(d$species_group, levels = species_levels)

theme_qc <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.background = element_rect(fill = "#F2F2F2", colour = "#B3B3B3"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, colour = "#4D4D4D")
  )

draw_object <- function(obj) {
  if (inherits(obj, "ggplot")) print(obj) else grid::grid.draw(obj)
}

save_all <- function(obj, stem, width, height) {
  png(file.path(output_dir, paste0(stem, ".png")), width = width * 300,
      height = height * 300, res = 300)
  draw_object(obj); dev.off()
  pdf(file.path(output_dir, paste0(stem, ".pdf")), width = width, height = height,
      useDingbats = FALSE)
  draw_object(obj); dev.off()
  svg(file.path(output_dir, paste0(stem, ".svg")), width = width, height = height)
  draw_object(obj); dev.off()
}

# Figure 1: BUSCO composition and internal-stop fraction.
n <- nrow(d)
busco_long <- data.frame(
  display_label = rep(d$display_label, times = 4),
  species_group = rep(d$species_group, times = 4),
  category = factor(rep(c("Single-copy", "Duplicated", "Fragmented", "Missing"), each = n),
                    levels = c("Missing", "Fragmented", "Duplicated", "Single-copy")),
  percent = c(d$busco_single_percent, d$busco_duplicated_percent,
              d$busco_fragmented_percent, d$busco_missing_percent)
)
busco_cols <- c("Single-copy" = "#0072B2", "Duplicated" = "#56B4E9",
                "Fragmented" = "#E69F00", "Missing" = "#D55E00")

p1a <- ggplot(busco_long, aes(percent, display_label, fill = category)) +
  geom_col(width = 0.82) +
  facet_grid(species_group ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = busco_cols) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 100, 20)) +
  coord_cartesian(xlim = c(0, 100)) +
  labs(title = "BUSCO completeness across 28 Aspergillus assemblies",
       subtitle = "Eurotiales odb12.2; genome mode with Miniprot",
       x = "BUSCO groups (%)", y = NULL, fill = NULL) + theme_qc

p1b <- ggplot(d, aes(busco_internal_stop_percent, display_label, colour = species_group)) +
  geom_point(size = 2.2, alpha = 0.9) +
  facet_grid(species_group ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = palette_values, guide = "none") +
  labs(title = "Complete BUSCOs with internal stops", x = "Internal-stop BUSCOs (%)", y = NULL) +
  theme_qc + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

fig1 <- arrangeGrob(p1a, p1b, ncol = 2, widths = c(3.2, 1.2))
save_all(fig1, "figure1_busco_completeness", 13, 12)

# Reusable horizontal point-plot helper.
metric_plot <- function(column, title, xlab, transform = "identity", formatter = waiver()) {
  ggplot(d, aes(x = .data[[column]], y = display_label, colour = species_group)) +
    geom_point(size = 2.1, alpha = 0.9) +
    scale_colour_manual(values = palette_values, guide = "none") +
    scale_x_continuous(trans = transform, labels = formatter) +
    labs(title = title, x = xlab, y = NULL) + theme_qc
}

# Figure 2: assembly composition and contiguity.
p2a <- metric_plot("total_length_bp", "Assembly size", "Genome size (Mb)", formatter = label_number(scale = 1e-6, accuracy = 0.1))
p2b <- metric_plot("sequences", "Scaffold count", "Scaffolds (log scale)", transform = "log10", formatter = label_comma()) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
p2c <- metric_plot("N50_bp", "Scaffold N50", "N50 (kb)", transform = "log10", formatter = label_number(scale = 1e-3, accuracy = 1))
p2d <- metric_plot("GC_percent", "GC composition", "GC (%)", formatter = label_number(accuracy = 0.1)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
fig2 <- arrangeGrob(p2a, p2b, p2c, p2d, ncol = 2,
                    top = grid::textGrob("Assembly size, contiguity and composition", gp = grid::gpar(fontface = "bold", fontsize = 14)))
save_all(fig2, "figure2_assembly_metrics", 13, 13)

# Figure 3: read support.
p3a <- metric_plot("mapped_percent", "Reads mapped", "Mapped reads (%)", formatter = label_number(accuracy = 0.1))
p3b <- metric_plot("properly_paired_percent", "Properly paired", "Properly paired reads (%)", formatter = label_number(accuracy = 0.1)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
p3c <- metric_plot("coverage_breadth_percent", "Assembly breadth covered", "Covered bases (%)", formatter = label_number(accuracy = 0.01))
p3d <- metric_plot("mean_depth", "Mean sequencing depth", "Mean depth (×)", formatter = label_number(accuracy = 1)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
fig3 <- arrangeGrob(p3a, p3b, p3c, p3d, ncol = 2,
                    top = grid::textGrob("Read support for each assembly", gp = grid::gpar(fontface = "bold", fontsize = 14)))
save_all(fig3, "figure3_read_support", 13, 13)

# Figure 4: integrated relationships with objective outlier labels.
d$outlier_label <- ifelse(d$busco_complete_percent < 95, as.character(d$sample), "")
p4a <- ggplot(d, aes(N50_bp / 1e3, busco_complete_percent, colour = species_group, size = total_length_bp / 1e6)) +
  geom_point(alpha = 0.85) +
  geom_text(aes(label = outlier_label), nudge_x = 0.5, nudge_y = 1.5, hjust = 0,
            size = 3, check_overlap = TRUE, show.legend = FALSE) +
  scale_colour_manual(values = palette_values) +
  scale_x_log10(labels = label_number(accuracy = 1)) +
  scale_size_continuous(range = c(2, 6)) +
  guides(size = "none") +
  labs(title = "Completeness versus contiguity", x = "Scaffold N50 (kb; log scale)",
       y = "Complete BUSCOs (%)", colour = "Species group", size = "Genome size (Mb)") + theme_qc

p4b <- ggplot(d, aes(mapped_percent, busco_complete_percent, colour = species_group, size = mean_depth)) +
  geom_point(alpha = 0.85) +
  geom_text(aes(label = outlier_label), nudge_x = 0.05, nudge_y = 1.5, hjust = 0,
            size = 3, check_overlap = TRUE, show.legend = FALSE) +
  scale_colour_manual(values = palette_values) +
  scale_size_continuous(range = c(2, 6)) +
  guides(size = "none", colour = "none") +
  labs(title = "Completeness versus read support", x = "Mapped reads (%)",
       y = "Complete BUSCOs (%)", colour = "Species group", size = "Mean depth (×)") + theme_qc

fig4 <- arrangeGrob(p4a, p4b, ncol = 2,
                    top = grid::textGrob("Integrated assembly QC overview", gp = grid::gpar(fontface = "bold", fontsize = 14)))
save_all(fig4, "figure4_integrated_qc", 13, 6.5)

# Machine-readable figure data and concise notes.
write.csv(d, file.path(output_dir, "figure_source_data.csv"), row.names = FALSE, na = "")
notes <- c(
  "# Aspergillus assembly QC figure drafts",
  "",
  "These figures summarize the 28 final SPAdes `scaffolds.fasta` assemblies.",
  "",
  "- Figure 1: BUSCO composition and the fraction of complete BUSCO predictions containing internal stop codons.",
  "- Figure 2: genome size, scaffold count, N50 and GC content.",
  "- Figure 3: mapped reads, proper pairs, breadth and mean depth.",
  "- Figure 4: integrated completeness/contiguity and completeness/read-support views.",
  "",
  "Assemblies below 95% complete BUSCO are labelled in Figure 4 for visibility; this is an annotation convention, not an automatic exclusion rule.",
  "PNG files are suitable for quick review; PDF and SVG are editable vector versions.",
  "The plotting script and source CSV are included for reproducibility."
)
writeLines(notes, file.path(output_dir, "README.md"))
