Human_Mouse_Counts_Analysis
================

# Here we are importing data on human and mouse ortholog genes:

``` r
data <- read.csv("/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/genehomies/RESULTS/HUMAN/human_gene_names.csv")

sig_genes <- read.csv("/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/genehomies/RESULTS/HUMAN/sig_genes_human.csv")

# significant ortholog human/mouse
load("/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/genehomies/RESULTS/MOUSEVSHUMAN/orthologs_results.RData")
```

# Creating dataframes we will use to plot counts data in human & mouse.

``` r
#Load counts data for mouse
MS_counts <- read.table("/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/genehomies/DATA/MOUSE/salmon.merged.gene_counts.tsv", header=TRUE, row.names=1)
#Load counts data for human
human_counts <- read.table("/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/genehomies/DATA/HUMAN/salmon.merged.gene_counts.tsv", header=TRUE, row.names=1)

human_brain_counts <- human_counts %>%
  filter(gene_name %in% list("HTR5A", "IER3"))

genes_of_interest <- conserved_changing_orthologues_dereplicated

human_genes_of_interest_counts <- human_counts %>%
  filter(gene_name %in% genes_of_interest$Gene.name)

human_gene_data_list <- lapply(human_brain_counts$gene_name, function(gene) {
  human_brain_counts[human_brain_counts$gene_name == gene, ]
})

human_split <- strsplit(colnames(human_gene_data_list[[1]]), '_')


mouse_brain_counts <- MS_counts %>%
  filter(gene_name %in% list("Htr5a", "Ier3"))

mouse_genes_of_interest_counts <- MS_counts %>%
  filter(gene_name %in% genes_of_interest$Mouse.gene.name)

mouse_gene_data_list <- lapply(mouse_brain_counts$gene_name, function(gene) {
  mouse_brain_counts[mouse_brain_counts$gene_name == gene, ]
})


# Define your time points as strings (to match column names exactly)
time_points <- c("0", "0.5", "1", "2", "2.5", "3", "3.5", "4", "4.5", "5", "5.5", "12", "24", "48", "96")

# Pick one gene for now (apply this to all later)
gene_df <- human_gene_data_list[[1]]

# Initialize list to store avg/sd for each timepoint
human_avg_sd_list <- list()

# Loop through each time point and calculate stats
for (tp in time_points) {
  # Match columns like gfp_<tp>_rep
  cols <- grep(paste0("gfp_", tp, "_"), colnames(gene_df))
  
  if (length(cols) > 0) {
    avg <- rowMeans(gene_df[, cols])
    std_dev <- apply(gene_df[, cols], 1, sd)
    human_avg_sd_list[[tp]] <- data.frame(avg, std_dev)
  }
}

# Combine all into one dataframe
combined_avg_sd <- do.call(cbind, human_avg_sd_list)

# Rename columns: timepoint_avg, timepoint_sd
colnames(combined_avg_sd) <- as.vector(t(outer(names(human_avg_sd_list), c("avg", "sd"), paste, sep = "_")))

# View result
head(combined_avg_sd)
```

    ##                       0_avg     0_sd 0.5_avg 0.5_sd    1_avg     1_sd    2_avg
    ## ENSG00000137331.12 804.8152 267.8713 856.267  328.7 770.7345 142.6765 754.0355
    ##                        2_sd  2.5_avg   2.5_sd    3_avg     3_sd 3.5_avg
    ## ENSG00000137331.12 41.72142 750.7325 72.99958 705.0725 14.58125  639.92
    ##                      3.5_sd    4_avg     4_sd  4.5_avg   4.5_sd    5_avg
    ## ENSG00000137331.12 99.03738 563.7255 95.81226 710.5485 14.11032 490.4405
    ##                        5_sd  5.5_avg   5.5_sd   12_avg    12_sd   24_avg
    ## ENSG00000137331.12 51.74678 568.1245 12.41609 235.8963 56.88562 313.0523
    ##                       24_sd   48_avg    48_sd   96_avg    96_sd
    ## ENSG00000137331.12 36.16178 217.0053 52.41033 289.4453 46.68451

# Creating split x-axis plot for viewing counts expression across all timepoints in human brain genes.

Genes of interest: Htr5a, Ier3 Timepoints: 0, 0.5, 1, 1.5, 2, 2.5, 3,
3.5, 4, 4.5, 5, 5.5, 12, 24, 48, 92

``` r
# Define the split point
split_point <- 5.5

# Prepare data for plotting
human_df_list <- list()
for (gene in unique(human_brain_counts$gene_name)) {
  gene_df <- human_brain_counts %>% filter(gene_name == gene)
  human_avg_sd_list <- list()
  for (tp in time_points) {
    cols <- grep(paste0("gfp_", tp, "_"), colnames(gene_df))
    if (length(cols) > 0) {
      avg <- rowMeans(gene_df[, cols])
      std_dev <- apply(gene_df[, cols], 1, sd)
      human_avg_sd_list[[tp]] <- data.frame(timepoint = tp, mean = avg, sd = std_dev)
    }
  }
  if (length(human_avg_sd_list) > 0) {
    combined_gene_avg_sd <- bind_rows(human_avg_sd_list)
    human_df_list[[gene]] <- combined_gene_avg_sd
  }
}

# Loop through each gene and create split axis plots
for (gene in names(human_df_list)) {
  df <- human_df_list[[gene]]
  df$timepoint <- as.numeric(as.character(df$timepoint))

  # Filter data for the two sections of the plot
  df_early <- df %>% filter(timepoint <= split_point)
  df_late <- df %>% filter(timepoint >= split_point)
  
  y_min <- min(c(df_early$mean - df_early$sd, df_late$mean - df_late$sd))
  y_max <- max(c(df_early$mean + df_early$sd, df_late$mean + df_late$sd))

  # Create the first plot (early time points)
  p1 <- ggplot(df_early, aes(x = timepoint, y = mean, group = 1)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
    scale_x_continuous(breaks = unique(df_early$timepoint), limits = c(min(df_early$timepoint), split_point)) +
    scale_y_continuous(limits = c(y_min, y_max)) +
    labs(
      title = paste(gene, "Expression Across Time (Human)"),
      y = "Mean Count",
      x = "Time (hours)"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title.x = element_blank() # Remove x-axis title for the first plot
    )

  # Create the second plot (later time points)
  p2 <- ggplot(df_late, aes(x = timepoint, y = mean, group = 1)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
    scale_x_continuous(breaks = unique(df_late$timepoint), limits = c(split_point, max(df_late$timepoint))) + # Corrected line
    scale_y_continuous(limits = c(y_min, y_max)) +
    labs(
      y = NULL, # Remove y-axis label for the second plot
      x = "Time (hours)"
    ) +
    theme(
      axis.title.y = element_blank(), # Remove y-axis title for the second plot
      axis.text.y = element_blank()    # Remove y-axis text for the second plot
    )

  # Combine the two plots using patchwork
  combined_plot <- p1 + p2 + plot_layout(ncol = 2, widths = c(2, 3)) # Adjust widths as needed

  # Add a common caption
  combined_plot <- combined_plot + plot_annotation(caption = "Error bars represent standard deviation") &
    theme(plot.caption = element_text(hjust = 0.5))

  # Save the combined plot
  ggsave(filename = paste0("figures/", gene, "_human_split_expression.png"), plot = combined_plot, width = 8, height = 4)
  print(combined_plot)
}
```

![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20human%20brain%20genes-1.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20human%20brain%20genes-2.png)<!-- -->

``` r
# Now for the combined facet plot with a split axis (more complex)
# We'll need to manipulate the data and plotting for this

# Create a function to generate the split facet plot for a single gene
split_facet_plot_gene <- function(gene_data) {
  split_point <- 5.5
  gene_name <- unique(gene_data$gene)

  df_early <- gene_data %>% filter(timepoint <= split_point)
  df_late <- gene_data %>% filter(timepoint >= split_point)
  
  y_min <- min(c(df_early$mean - df_early$sd, df_late$mean - df_late$sd))
  y_max <- max(c(df_early$mean + df_early$sd, df_late$mean + df_late$sd))


  p1 <- ggplot(df_early, aes(x = timepoint, y = mean, group = 1)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
    scale_x_continuous(breaks = unique(df_early$timepoint), limits = c(min(df_early$timepoint), split_point)) +
    scale_y_continuous(limits = c(y_min, y_max)) + 
    labs(y = "Mean Count") +
    theme(axis.title.x = element_blank())

  p2 <- ggplot(df_late, aes(x = timepoint, y = mean, group = 1)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
    scale_x_continuous(breaks = unique(df_late$timepoint), limits = c(split_point, max(df_late$timepoint))) + # Corrected line
    scale_y_continuous(limits = c(y_min, y_max)) + 
    theme(axis.title.y = element_blank(), axis.text.y = element_blank())

  combined <- p1 + p2 + plot_layout(ncol = 2, widths = c(2, 3)) +
    plot_annotation(title = gene_name) &
    theme(plot.title = element_text(hjust = 0.5))

  return(combined)
}
```

``` r
# Define the split point
split_point <- 5.5

# Prepare data for plotting
all_human_df_list <- list()
for (gene in unique(human_genes_of_interest_counts$gene_name)) {
  all_gene_df <- human_genes_of_interest_counts %>% filter(gene_name == gene)
  all_human_avg_sd_list <- list()
  for (tp in time_points) {
    cols <- grep(paste0("gfp_", tp, "_"), colnames(all_gene_df))
    if (length(cols) > 0) {
      avg <- rowMeans(all_gene_df[, cols])
      std_dev <- apply(all_gene_df[, cols], 1, sd)
      all_human_avg_sd_list[[tp]] <- data.frame(timepoint = tp, mean = avg, sd = std_dev)
    }
  }
  if (length(human_avg_sd_list) > 0) {
    combined_gene_avg_sd <- bind_rows(all_human_avg_sd_list)
    all_human_df_list[[gene]] <- combined_gene_avg_sd
  }
}

# Loop through each gene and create split axis plots
for (gene in names(all_human_df_list)) {
  df <- all_human_df_list[[gene]]
  df$timepoint <- as.numeric(as.character(df$timepoint))

  # Filter data for the two sections of the plot
  df_early <- df %>% filter(timepoint <= split_point)
  df_late <- df %>% filter(timepoint >= split_point)
  
  y_min <- min(c(df_early$mean - df_early$sd, df_late$mean - df_late$sd))
  y_max <- max(c(df_early$mean + df_early$sd, df_late$mean + df_late$sd))

  # Create the first plot (early time points)
  p1 <- ggplot(df_early, aes(x = timepoint, y = mean, group = 1)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
    scale_x_continuous(breaks = unique(df_early$timepoint), limits = c(min(df_early$timepoint), split_point)) +
    scale_y_continuous(limits = c(y_min, y_max)) +
    labs(
      title = paste("    ", gene, "Expression Across Time (Human)"),
      y = "Mean Count",
      x = "Time (hours)"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title.x = element_blank() # Remove x-axis title for the first plot
    )

  # Create the second plot (later time points)
  p2 <- ggplot(df_late, aes(x = timepoint, y = mean, group = 1)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
    scale_x_continuous(breaks = unique(df_late$timepoint), limits = c(split_point, max(df_late$timepoint))) + # Corrected line
    scale_y_continuous(limits = c(y_min, y_max)) +
    labs(
      y = NULL, # Remove y-axis label for the second plot
      x = "Time (hours)"
    ) +
    theme(
      axis.title.y = element_blank(), # Remove y-axis title for the second plot
      axis.text.y = element_blank()    # Remove y-axis text for the second plot
    )

  # Combine the two plots using patchwork
  combined_plot <- p1 + p2 + plot_layout(ncol = 2, widths = c(2, 3)) # Adjust widths as needed

  # Add a common caption
  combined_plot <- combined_plot + plot_annotation(caption = "Error bars represent standard deviation") &
    theme(plot.caption = element_text(hjust = 0.5))

  # Save the combined plot
  ggsave(filename = paste0("figures/", gene, "_human_split_expression.png"), plot = combined_plot, width = 8, height = 4)
  print(combined_plot)
}
```

![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-1.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-2.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-3.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-4.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-5.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-6.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-7.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-8.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-9.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-10.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-11.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-12.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-13.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-14.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-15.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-16.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-17.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-18.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-19.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-20.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-21.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-22.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-23.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-24.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-25.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-26.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-27.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-28.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-29.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-30.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-31.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-32.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-33.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-34.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-35.png)<!-- -->![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-36.png)<!-- -->

``` r
# Now for the combined facet plot with a split axis (more complex)
# We'll need to manipulate the data and plotting for this

# Combine all gene data for a facet plot
all_human_summary <- dplyr::bind_rows(all_human_df_list, .id = "gene")
all_human_summary$timepoint <- as.numeric(as.character(all_human_summary$timepoint))

reduced_time_points <- c(0, 12, 24, 48, 96)
all_human_summary <- all_human_summary[all_human_summary$timepoint %in% reduced_time_points, ]

# Create a function to generate the split facet plot for a single gene
facet_plot <- ggplot(all_human_summary, aes(x = as.numeric(timepoint), y = mean, group = gene)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
  scale_x_continuous(breaks = unique(all_human_summary$timepoint)) +
  facet_wrap(~ gene, scales = "free_y") +
  labs(
    title = "Expression Across Time(Human Counts)",
    y = "Mean Count",
    x = "Time (hours)",
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )

# save the facet plot to a file
ggsave(filename = "figures/all_genes_facet_expression.png",
       plot = facet_plot,
       width = 12,
       height = 8)
print(facet_plot)
```

![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-37.png)<!-- -->

``` r
# Prepare data for plotting
mouse_df_list <- list()
mouse_time_points <- c("0", "12", "24","48", "96")
for (gene in unique(mouse_brain_counts$gene_name)) {
  mouse_gene_df <- mouse_brain_counts %>% filter(gene_name == gene)
  mouse_avg_sd_list <- list()
  for (tp in mouse_time_points) {
    cols <- grep(paste0("WT_", tp, "_"), colnames(mouse_gene_df))
    if (length(cols) > 0) {
      avg <- rowMeans(gene_df[, cols])
      std_dev <- apply(gene_df[, cols], 1, sd)
      mouse_avg_sd_list[[tp]] <- data.frame(timepoint = tp, mean = avg, sd = std_dev)
    }
  }
  if (length(mouse_avg_sd_list) > 0) {
    combined_gene_avg_sd <- bind_rows(mouse_avg_sd_list)
    mouse_df_list[[gene]] <- combined_gene_avg_sd
  }
}

for (gene in names(mouse_df_list)) {
  df <- mouse_df_list[[gene]]
  df$timepoint <- as.numeric(as.character(df$timepoint)) 
  p <- ggplot(df, aes(x = timepoint, y = mean, group = 1)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
    scale_x_continuous(breaks = unique(df$timepoint)) +
    labs(
      title = paste(gene, "Mouse Expression Across Time"),
      y = "Mean Count",
      x = "Time (hours)",
      caption = "Error bars represent standard error of the mean"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(face = "bold")
    )

# save the image in the figures folder
ggsave(filename = paste0("figures/", gene, "_mouse_expression.png"), plot = p, width = 6, height = 4)
}
```

``` r
# Define the split point
split_point <- 5.5

# Prepare data for plotting
all_mouse_df_list <- list()
for (gene in unique(mouse_genes_of_interest_counts$gene_name)) {
  all_mouse_gene_df <- mouse_genes_of_interest_counts %>% filter(gene_name == gene)
  all_mouse_avg_sd_list <- list()
  for (tp in mouse_time_points) {
    cols <- grep(paste0("WT_", tp, "_"), colnames(all_mouse_gene_df))
    if (length(cols) > 0) {
      avg <- rowMeans(all_mouse_gene_df[, cols])
      std_dev <- apply(all_mouse_gene_df[, cols], 1, sd)
      all_mouse_avg_sd_list[[tp]] <- data.frame(timepoint = tp, mean = avg, sd = std_dev)
    }
  }
  if (length(all_mouse_avg_sd_list) > 0) {
    mouse_combined_gene_avg_sd <- bind_rows(all_mouse_avg_sd_list)
    all_mouse_df_list[[gene]] <- mouse_combined_gene_avg_sd
  }
}

# Loop through each gene and create split axis plots
for (gene in names(all_mouse_df_list)) {
  df <- all_mouse_df_list[[gene]]
  df$timepoint <- as.numeric(as.character(df$timepoint))

  p <- ggplot(df, aes(x = timepoint, y = mean, group = 1)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
    scale_x_continuous(breaks = unique(df$timepoint)) +
    labs(
      title = paste(gene, "Mouse Expression Across Time"),
      y = "Mean Count",
      x = "Time (hours)",
      caption = "Error bars represent standard error of the mean"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(face = "bold")
    )


  # Save the combined plot
  ggsave(filename = paste0("figures/", gene, "_mouse_expression.png"), plot = p, width = 8, height = 4)
}

# Now for the combined facet plot with a split axis (more complex)
# We'll need to manipulate the data and plotting for this

# Combine all gene data for a facet plot
all_mouse_summary <- dplyr::bind_rows(all_mouse_df_list, .id = "gene")
all_mouse_summary$timepoint <- as.numeric(as.character(all_mouse_summary$timepoint))

# Create a function to generate the split facet plot for a single gene
facet_plot <- ggplot(all_mouse_summary, aes(x = as.numeric(timepoint), y = mean, group = gene)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
  scale_x_continuous(breaks = unique(all_mouse_summary$timepoint)) +
  facet_wrap(~ gene, scales = "free_y") +
  labs(
    title = "Expression Across Time(Mouse Counts)",
    y = "Mean Count",
    x = "Time (hours)",
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )

# save the facet plot to a file
ggsave(filename = "figures/all_mouse_counts_facet_expression.png",
       plot = facet_plot,
       width = 12,
       height = 8)
print(facet_plot)
```

![](human_mouse_analysis_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20mouse%20genes-1.png)<!-- -->
