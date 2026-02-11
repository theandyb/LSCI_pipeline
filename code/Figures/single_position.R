library(tidyverse)
library(ggpubfigs)
library(ggsci)

source("code/Figures/load_data.R")


# This loads the results from single position models for all subtypes for the ALL population
df_all <- load_all_pos_results(c(subtypes, "all_GC_AT"))
df_all$cramer <- sqrt(df_all$dev / (df_all$singletons + df_all$controls))
df_all$side <- ifelse(df_all$offset > 0, 1, -1)

# Plot log relative entropys at all positions for each subtype
df_all %>%
    filter(subtype == "AT_GC") %>%
    mutate(is_sig = pchisq(dev, 3, lower.tail = FALSE) < 0.05) %>% # change order to TRUE then FALSE
    mutate(is_sig = factor(is_sig, levels = c(TRUE, FALSE), labels = c("True", "False"))) %>%
    ggplot(aes(x = offset, y = log10(re), color = is_sig)) +
    geom_point() +
    theme_classic() +
    scale_color_npg() + 
    xlab("Flanking Position") +
    ylab("Log10 Relative Entropy") +
    labs(color = "Significant\n(p < 0.05)") + # bigger axis text
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 22),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 18),
          legend.position = c(0.8, 0.8)) +
    guides(color = guide_legend(override.aes = list(shape = 15, size = 3)))

# Plot overall statistics and residuals
resid_logo_plot <- function(res_dir, pop, st, df_all, rp = c(-5:-1,1:5)){
  res_df <- get_pos_resid_all(res_dir, pop, st, rp)

  mat_df <- res_df %>%
    filter(status == "singletons") %>%
    rowwise() %>%
    mutate(re = sign(res)*re_res) %>%
    select(rp, nuc, re) %>%
    pivot_wider(names_from = nuc, values_from = re) %>%
    bind_rows(data.frame(rp=0, A=0, C=0, G=0, T=0)) %>%
    arrange(rp) %>%
    select(-rp) %>%
    as.data.frame()

  p <- ggseqlogo::ggseqlogo(t(mat_df), method = "custom") +
  theme_classic() +
    scale_x_continuous(breaks = 1:11, labels = function(x){x - 6}) +
    geom_hline(yintercept = 0, linetype = "dotted")
  return(p)
}

resid_combo_plot <- function(res_dir, pop, st, df_all, rp = c(-5:-1,1:5)){
  res_df <- get_pos_resid_all(res_dir, pop, st, rp)
  res_df$ID <- as.character(res_df$rp)
  xlab <- as.character(min(res_df$rp):max(res_df$rp))
  res_df$ID <- factor(res_df$ID, levels = xlab)

  df_stat <- df_all %>%
    filter(subtype == st) %>%
    arrange(offset) %>%
    filter(offset %in% rp)

  df_stat$ID <- factor(df_stat$offset, levels = xlab)
  df_stat$group <- ifelse(df_stat$offset < 0, -1, 1)

  res_df2 <- res_df %>%
    mutate(sign_re = sign(res) * re_res) %>%
    select(ID, nuc, status, sign_re) %>%
    pivot_wider(names_from = status, values_from = sign_re) %>%
    mutate(re_cont = sign(singletons) * (abs(singletons) + abs(controls))) 

  p <- res_df2 %>%
    ggplot(aes(x= ID, y = re_cont, fill = nuc)) +
    scale_x_discrete(labels = xlab, drop = FALSE) +
    geom_bar(position = "stack", stat = "identity") +
    geom_point(aes(y = re, fill = NULL), data = df_stat, show.legend = F, size = 3) +
    geom_line(aes(y = re, fill = NULL, group = group), 
              linewidth = 1, data = df_stat, linetype = "dotted", alpha = 1, colour = "grey") +
    scale_fill_manual(values =ggpubfigs::friendly_pals$nickel_five) +
    ylab("") +
    theme_classic() +
    xlab(paste0(subtype_print_names(st))) +
    labs(fill = "") +
    theme(legend.position = "none") + # make axis text bigger
    theme(axis.text = element_text(size = 18),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 22, hjust = 0.5))
  return(p)
}

resid_combo_plot(single_pos_dir, "ALL", "AT_CG", df_all, rp = c(-5:-1,1:5))

# Get all subtypes' plots and store in a list
plot_list <- vector(mode = "list", length = 9)
names(plot_list) <- subtypes

for(st in subtypes){
  plot_list[[st]] <- resid_combo_plot(single_pos_dir, "ALL", st, df_all, rp = c(-5:-1,1:5))
}

plot_list[["AT_CG"]] <- plot_list[["AT_CG"]] +
  theme(legend.position = "inside",
          legend.position.inside = c(0.2, 0.9), 
          legend.direction = "horizontal",
          legend.text = element_text(size = 8), 
          legend.key.size = unit(14, "pt")) 

# 1. Find the maximum y-axis limit across all plots
max_y <- -Inf
min_y <- Inf
for (plot in plot_list) {
  y_range <- layer_scales(plot)$y$range$range
  if (!is.null(y_range) && y_range[2] > max_y) {
    max_y <- y_range[2]
  }
  if (!is.null(y_range) && y_range[1] < min_y) {
    min_y <- y_range[1]
  }
}

# 2. Update all plots in the list to have the same y-axis limits
updated_plot_list <- lapply(plot_list, function(plot) {
  plot + ylim(min_y - 0.001, max_y + 0.001) # Assuming you want the y-axis to start at 0
})

#sjPlot::plot_grid(x = updated_plot_list, margin = c(0.1,0.1,0.1,0.1))
gridExtra::grid.arrange(grobs = updated_plot_list, ncol = 3)
