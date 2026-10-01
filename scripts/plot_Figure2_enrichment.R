library(dplyr)
library(ggplot2)
library(readr)


## Input files only; no randomization is performed here.
summary_constrained <- read.delim("bootstrap_constrained_trans_eqtl_summary_1000.txt")
hic_summary <- read_tsv("hic_constrained_random_summary.tsv", show_col_types = FALSE)

category_cols <- c(
  "Inter-chromosomal\nnon-syntenic" = "#F3B5B3",
  "Inter-chromosomal\nsyntenic" = "#F5E69A",
  "Intra-chromosome\ndifferent syntenic block" = "#B5AFDE",
  "Intra-chromosome\nsame syntenic block" = "#B6D8B6"
)


plot_df <- summary_constrained %>%
  mutate(
    category_label = recode(
      duplicatestatus,
      "Inter-chromosomal non-syntenic" = "Inter-chromosomal\nnon-syntenic",
      "Inter-chromosomal syntenic" = "Inter-chromosomal\nsyntenic",
      "Intra-chromosome different syntenic block" = "Intra-chromosome\ndifferent syntenic block",
      "Intra-chromosome same syntenic block" = "Intra-chromosome\nsame syntenic block"
    ),
    category_label = factor(category_label, levels = names(category_cols)[c(4, 3, 2, 1)]),
    label_p = ifelse(mean_enrichment >= 1, empirical_p_upper, empirical_p_lower),
    p_label = paste0("p=", format.pval(label_p, digits = 2, eps = 0.001))
  )


p2d <- ggplot(
  plot_df,
  aes(category_label, mean_enrichment, fill = category_label)
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_col(width = 0.65) +
  geom_errorbar(
    aes(
      ymin = mean_enrichment - se_enrichment,
      ymax = mean_enrichment + se_enrichment
    ),
    width = 0.18
  ) +
  geom_text(
    aes(
      label = observed_n,
      y = mean_enrichment + 0.025
    ),
    vjust = 0,
    size = 5.8
  ) +
  geom_text(
    aes(
      label = p_label,
      y = mean_enrichment + 0.18
    ),
    vjust = 0,
    size = 4.8
  ) +
  scale_fill_manual(
    values = category_cols,
    guide = "none"
  ) +
  scale_y_continuous(
    limits = c(0, 1.65),
    breaks = seq(0, 1.6, 0.2),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = NULL,
    y = "Enrichment over constrained randomization"
  ) +
  theme_classic(base_size = 20) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 14),
    plot.margin = margin(8, 8, 8, 8)
  )


p2d


ggsave("Figure2D_constrained_enrichment_horizontal.png", p2d, width = 7.2, height = 3.6, dpi = 600)
ggsave("Figure2D_constrained_enrichment_horizontal.pdf", p2d, width = 7.2, height = 3.6)


hic_plot_df <- hic_summary %>%
  mutate(
    category_label = recode(
      Category,
      "Inter-chromosomal non-syntenic" = "Inter-chromosomal\nnon-syntenic",
      "Inter-chromosomal syntenic" = "Inter-chromosomal\nsyntenic",
      "Intra-chromosome different syntenic block" = "Intra-chromosome\ndifferent syntenic block",
      "Intra-chromosome same syntenic block" = "Intra-chromosome\nsame syntenic block"
    ),
    category_label = factor(category_label, levels = names(category_cols)[c(4, 3, 2, 1)]),
    p_label = empirical_p_contact_high,
    ratio_label = sprintf("%.2f", mean_contact_ratio),
    n_label = paste0("n=", obs_n),
    p_text = paste0("p=", format.pval(p_label, digits = 2, eps = 0.001))
  )


xmax <- max(hic_plot_df$mean_contact_ratio, na.rm = TRUE) * 1.35
p2e <- ggplot(
  hic_plot_df,
  aes(category_label, mean_contact_ratio, fill = category_label)
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_col(width = 0.65) +
  geom_text(
    aes(label = n_label),
    y = 0.05,
    vjust = 0,
    size = 3.8
  ) +
  geom_text(
    aes(label = ratio_label),
    nudge_y = 0.04,
    vjust = 0,
    size = 5.8
  ) +
  geom_text(
    aes(label = p_text),
    nudge_y = 0.28,
    vjust = 0,
    size = 4.8
  ) +
  scale_fill_manual(
    values = category_cols,
    guide = "none"
  ) +
  scale_y_continuous(
    limits = c(0, xmax),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    x = NULL,
    y = "Observed / randomized Hi-C contact"
  ) +
  theme_classic(base_size = 20) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 14),
    plot.margin = margin(8, 8, 8, 8)
  )


p2e

ggsave("Figure2_HiC_contact_ratio_horizontal.png", p2e, width = 8.2, height = 3.8, dpi = 600, limitsize = FALSE)
ggsave("Figure2_HiC_contact_ratio_horizontal.pdf", p2e, width = 8.2, height = 3.8, limitsize = FALSE)


