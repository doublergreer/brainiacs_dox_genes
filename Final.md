Analysis_of_Human_and_Mouse_Gene_Expressions
================
M. Cheishvili, M. Todd, T. Whittaker, R. Greer
2025-25-04

Our objective in this project is to compare gene expression between
human and mouse genes after doxycycline (dox) treatment. We used raw
counts and transcripts per million (TPM) data from RNA-seq datasets for
both species.

To begin, we imported the large data sets and used filtering techniques
to reduced our focus to the genes that changed significantly, based on
the counts and TPM data. Once the significant orthologous genes were
determined, we analyzed them using the mean and standard deviation
expression values at a variety of time points. Furthermore, we
determined whether the regulation of the sig genes using directionality
vectors on both raw data and log fold change, to ensure resiliency.

To analyze patterns in the changes, we utilized many figures and
Integrative Genomics Viewer (IGV). For the sig genes we found, we
utilized the KEGG database to import pathway information. We found that
many of the genes involved transduction signaling, metabolic pathways,
and disease-associated mechanisms.

Overall, our goal was to conduct a comprehensive comparative analysis of
human and mouse genes after doxycycline treatment.

Human versus Mouse Analysis of Counts

# DNA Counts Analysis:

## Importing data on human and mouse ortholog genes

``` r
data <- read.csv("/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/genehomies/RESULTS/HUMAN/human_gene_names.csv")

sig_genes <- read.csv("/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/genehomies/RESULTS/HUMAN/sig_genes_human.csv")

# significant ortholog human/mouse
load("/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/genehomies/RESULTS/MOUSEVSHUMAN/orthologs_results.RData")
```

## Creating dataframes we will use to plot counts data in human & mouse.

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
time_points <- c("0", "0.5", "1", "1.5", "2", "2.5", "3", "3.5", "4", "4.5", "5", "5.5", "12", "24", "48", "96")

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

    ##                       0_avg     0_sd 0.5_avg 0.5_sd    1_avg     1_sd  1.5_avg
    ## ENSG00000137331.12 804.8152 267.8713 856.267  328.7 770.7345 142.6765 785.7865
    ##                      1.5_sd    2_avg     2_sd  2.5_avg   2.5_sd    3_avg
    ## ENSG00000137331.12 75.42496 754.0355 41.72142 750.7325 72.99958 705.0725
    ##                        3_sd 3.5_avg   3.5_sd    4_avg     4_sd  4.5_avg
    ## ENSG00000137331.12 14.58125  639.92 99.03738 563.7255 95.81226 710.5485
    ##                      4.5_sd    5_avg     5_sd  5.5_avg   5.5_sd   12_avg
    ## ENSG00000137331.12 14.11032 490.4405 51.74678 568.1245 12.41609 235.8963
    ##                       12_sd   24_avg    24_sd   48_avg    48_sd   96_avg
    ## ENSG00000137331.12 56.88562 313.0523 36.16178 217.0053 52.41033 289.4453
    ##                       96_sd
    ## ENSG00000137331.12 46.68451

## Creating split x-axis plot for viewing counts expression across all timepoints in human brain genes.

Our human counts data sampled expression across a higher number of
timepoints than our mouse data, so we can use a split x-axis plot to
make all the smaller timepoint data visible while still showing
correctly scaled late timepoint data in the second half of the graph.

### Genes of interest: HTR5A, IER3

### Timepoints: 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 12, 24, 48, 96

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

![](Final_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20human%20brain%20genes-1.png)<!-- -->![](Final_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20human%20brain%20genes-2.png)<!-- -->

## Creating facet plot of all significant orthologues

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
}

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

![](Final_files/figure-gfm/creating%20split-x%20axis%20plot%20for%20all%20significant%20human%20genes-1.png)<!-- -->

## Creating plot for viewing counts expression in mouse brain genes.

### Genes of interest: Htr5a, Ier3 

### Timepoints: 0, 12, 24, 48, 96

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
print(p)
}
```

![](Final_files/figure-gfm/creating%20normal%20plot%20for%20mouse%20brain%20genes-1.png)<!-- -->![](Final_files/figure-gfm/creating%20normal%20plot%20for%20mouse%20brain%20genes-2.png)<!-- -->

## Human brain genes expression curves for comparison:

Despite the differences in x-axis densities, we can see that the human
and mouse brain genes counts follow the same trends when exposed to
doxycycline.
![humanhtr5a-compar](figures/IER3_human_split_expression.png)
![humanier3-compar](figures/HTR5A_human_split_expression.png)

``` r
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

# Loop through each gene and create plots
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

![](Final_files/figure-gfm/creating%20plot%20for%20all%20significant%20mouse%20genes-1.png)<!-- -->

## Human facet plot for comparison:

![humancounts-compar](figures/all_genes_facet_expression.png) Here we
can see that for most of the genes, the expression trends are similar
from human to mouse. A couple notable exceptions are CBR3 (downregulated
in human) and Cbr3 (upregulated in mouse) and DLL1 (upregulated in
human) and Dll1(downregulated in mouse). Other than these two, the 0 to
96 timepoint changes between human and mouse seem to be moving in the
same direction, though sometimes with different curve shapes.

Import the conserved human and mouse genes again if you cleared your
environment before starting this portion of the document.

## We will then load in the TPM datasets for both human and mouse

We then will proceed to filter both data sets by the 36 genes that were
found to be conserved in both human and mouse

In the human data there are additional time points and replicates, we
will remove those additional time points and columns and simply focus on
0, 12, 24, 48. 96 timepoints and keep the first 3 replicants

We will then calculate the average and standard deviations for both data
sets so that we can proceed to plot the data using ggplot

``` r
# Loading in the conserved gene data between  Human and Mouse 

data <- read.csv("/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/genehomies/RESULTS/HUMAN/human_gene_names.csv")

sig_genes <- read.csv("/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/genehomies/RESULTS/HUMAN/sig_genes_human.csv")

#significant ortholog human/mouse
load("/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/genehomies/RESULTS/MOUSEVSHUMAN/orthologs_results.RData")


#Load TPM data for mouse
MS_tpm <- read.table("/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/genehomies/DATA/MOUSE/salmon.merged.gene_tpm.tsv", header=TRUE, row.names=1)

#Load TPM data for human
human_tpm <- read.table("/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/genehomies/DATA/HUMAN/salmon.merged.gene_tpm.tsv", header=TRUE, row.names=1)

human_conserved_tpm <- human_tpm %>%
  filter(gene_name %in% list("HTR5A", "IER3","CCDC3","MEI1","RBM24","DPPA3","PSTPIP2","DLL1","BCHE","CPM","OTX2","SLC40A1","BIK","PTPRN","GCK","IER3","PPP1R3C","ELAVL4","OAS3","ABCC2","LHX5","PADI3","HTR5A","PSD4","GADL1","FGFR3","CYP1A1","CHAC1","ETV4","RHCG","SP5","ROBO4","ITGA10","LY6H","GAP43","IGF2","ATP1A4","CBR3"))


mouse_conserved_tpm <- MS_tpm %>%
  filter(gene_name %in% list("Htr5a", "Ier3","Cbr3","Ccdc3","Dppa3","Mei1","Rbm24","Pstpip2","Dll1","Bche","Cpm","Otx2","Slc40a1","Bik","Ptprn","Gck","Ppp1r3c","Elavl4","Oas3","Abcc2","Lhx5","Padi3","Psd4","Gadl1","Fgfr3","Cyp1a1","Chac1","Etv4","Rhcg","Sp5","Robo4","Itga10","Ly6h","Gap43","Igf2","Atp1a4"))

# Time and replicate values for RNAseq Data
time_points <- c("0", "12", "24", "48", "96")
replicates <- c("_1", "_2", "_3")

# initialize list for results

# Adding gene_id as a column header  for the meta data column

# Exact column names to retain in human_conserved_tpm
keep_columns <- c( 
 "gene_name", "gfp_0_1", "gfp_0_2", "gfp_0_3",
  "gfp_12_1", "gfp_12_2", "gfp_12_3",
  "gfp_24_1", "gfp_24_2", "gfp_24_3",
  "gfp_48_1", "gfp_48_2", "gfp_48_3",
  "gfp_96_1", "gfp_96_2", "gfp_96_3"
)

# Subset only those columns
filtered_human_tpm <- human_conserved_tpm[, keep_columns]

filtered_human_tpm <- cbind(gene_id = rownames(filtered_human_tpm), filtered_human_tpm)
rownames(filtered_human_tpm) <- NULL

mouse_conserved_tpm <-cbind(gene_id = rownames(mouse_conserved_tpm), mouse_conserved_tpm)
rownames(mouse_conserved_tpm) <- NULL

# Calculate average and standard deviation for MOUSE
average_and_stddev_values_MS_tpm <- list()
for (tp in time_points) {
  cols <- grep(paste0("WT_", tp, "_"), colnames(mouse_conserved_tpm))
  avg <- rowMeans(mouse_conserved_tpm[, cols])
  std_dev <- apply(mouse_conserved_tpm[, cols], 1, sd)
  std_dev <- data.frame(std_dev)
  combined_mouse <- cbind(avg, std_dev)
  average_and_stddev_values_MS_tpm <- c(average_and_stddev_values_MS_tpm, list(combined_mouse))
}
average_and_stddev_values_MS_tpm <- do.call(cbind, average_and_stddev_values_MS_tpm)
colnames(average_and_stddev_values_MS_tpm) <- paste0(rep(time_points, each = 2), c("_avg", "_sd"))

# Add gene_id back
average_and_stddev_values_MS_tpm <- cbind(
  gene_id = mouse_conserved_tpm$gene_id,
  gene_name = mouse_conserved_tpm$gene_name,
  average_and_stddev_values_MS_tpm
)

# Convert to long format for plotting
mouse_plot_df <- average_and_stddev_values_MS_tpm %>%
  pivot_longer(cols = -c(gene_id, gene_name), names_to = "Metric", values_to = "TPM") %>%
  separate(Metric, into = c("Time", "Type")) %>%
  pivot_wider(names_from = Type, values_from = TPM)

# Convert Time to a factor with correct order
mouse_plot_df$Time <- factor(mouse_plot_df$Time, levels = as.character(time_points))


# Add column names for the time points
colnames(average_and_stddev_values_MS_tpm) <- paste0(rep(time_points, each = 2), c("_avg", "_sd"))

# Repeat the above steps for the human tpm data 

# initialize list for results

# Calculate average and standard deviation for HUMAN
average_and_stddev_values_Human_tpm <- list()
for (tp in time_points) {
  cols <- grep(paste0("gfp_", tp, "_"), colnames(filtered_human_tpm))
  avg <- rowMeans(filtered_human_tpm[, cols])
  std_dev <- apply(filtered_human_tpm[, cols], 1, sd)
  std_dev <- data.frame(std_dev)
  combined_human <- cbind(avg, std_dev)
  average_and_stddev_values_Human_tpm <- c(average_and_stddev_values_Human_tpm, list(combined_human))
}
average_and_stddev_values_Human_tpm <- do.call(cbind, average_and_stddev_values_Human_tpm)
colnames(average_and_stddev_values_Human_tpm) <- paste0(rep(time_points, each = 2), c("_avg", "_sd"))

# Add gene_id back
average_and_stddev_values_Human_tpm <- cbind(
  gene_id = filtered_human_tpm$gene_id,
  gene_name = filtered_human_tpm$gene_name,
  average_and_stddev_values_Human_tpm
)

# Convert to long format for plotting
human_plot_df <- average_and_stddev_values_Human_tpm %>%
    pivot_longer(cols = -c(gene_id, gene_name), names_to = "Metric", values_to = "TPM") %>%
  separate(Metric, into = c("Time", "Type")) %>%
  pivot_wider(names_from = Type, values_from = TPM)

# Convert Time to a factor with correct order
human_plot_df$Time <- factor(human_plot_df$Time, levels = as.character(time_points))
```

``` r
# Mouse plot
mouse_viewable_plot <- ggplot(mouse_plot_df, aes(x = Time, y = avg, group = gene_id)) +
  geom_line() +
  geom_errorbar(aes(ymin = avg - sd, ymax = avg + sd), width = 0.2) +
  labs(x = 'Time (hours)', y = 'Average TPM') +
  ggtitle('Average TPM of Conserved Mouse Genes across Time Points') +
  facet_wrap(~ gene_name, scales = "free_y")
print(mouse_viewable_plot)
```

![](Final_files/figure-gfm/Plotting%20TPM%20Data%20from%20Mouse%20and%20Human%20Data-1.png)<!-- -->

``` r
# Human plot
human_viewable_plot <- ggplot(human_plot_df, aes(x = Time, y = avg, group = gene_id)) +
  geom_line() +
  geom_errorbar(aes(ymin = avg - sd, ymax = avg + sd), width = 0.2) +
  labs(x = 'Time (hours)', y = 'Average TPM') +
  ggtitle('Average TPM of Conserved Human Genes across Time Points') +
  facet_wrap(~ gene_name, scales = "free_y")

print(human_viewable_plot)
```

![](Final_files/figure-gfm/Plotting%20TPM%20Data%20from%20Mouse%20and%20Human%20Data-2.png)<!-- -->

From the line plots we see some similarities and differences between the
trends of these 36 genes in humans and mice. For the GCK, Cyp1a1 and
IER3 plots we see similar trends of downregulation in both humans and
mice. However in the majority of the genes we see there are differences
in the trends of upregulation and downregulation. For example, Cbr3 and
Oas3 are upregulated and downregulated respectively across the time
points. While these same two genes are downregulated and upregulated in
humans.

We will now use the TPM data to construct direction vectors for the
human and mouse datasets. This will allow us to see at what time points
the human and mouse gene expression transcripts match and how it
compares to the trends seen in the line plots illustrated above.

``` r
# Get intersecting genes using case-insensitive match
human_genes <- tolower(filtered_human_tpm$gene_name)
mouse_genes <- tolower(mouse_conserved_tpm$gene_name)
common_genes <- intersect(human_genes, mouse_genes)

# Ensure time points are clean integers
time_points <- c(0, 12, 24, 48, 96)


# Initialize result frame
n_transitions <- length(time_points) - 1
transition_labels <- paste0(
  "T", as.integer(time_points[-length(time_points)]),
  "_to_T", as.integer(time_points[-1])
)

direction_match <- data.frame(
  gene_name = common_genes,
  match_score = NA_integer_,
  total_transitions = n_transitions,
  matrix(NA_integer_, nrow = length(common_genes), ncol = n_transitions)
)

# Add trend columns OUTSIDE the loop
direction_match$human_trend <- NA_character_
direction_match$mouse_trend <- NA_character_

# Name time transition columns
colnames(direction_match)[4:ncol(direction_match)] <- transition_labels

# Loop over genes
for (i in seq_along(common_genes)) {
  gene <- common_genes[i]

  # Get expression rows
  h_row <- filtered_human_tpm[tolower(filtered_human_tpm$gene_name) == gene, ]
  m_row <- mouse_conserved_tpm[tolower(mouse_conserved_tpm$gene_name) == gene, ]

  # Extract average TPMs
  h_expr <- sapply(time_points, function(tp) {
    rowMeans(h_row[, grep(paste0("gfp_", tp, "_"), colnames(h_row)), drop = FALSE])
  })
  m_expr <- sapply(time_points, function(tp) {
    rowMeans(m_row[, grep(paste0("WT_", tp, "_"), colnames(m_row)), drop = FALSE])
  })

  # Get direction vectors
  h_dir <- sign(diff(h_expr))
  m_dir <- sign(diff(m_expr))

  # Compare direction
 match_vector <- character(length = n_transitions)
for (j in seq_len(n_transitions)) {
  if (h_dir[j] == m_dir[j]) {
    match_vector[j] <- "Match"
  } else if (h_dir[j] == 1 & m_dir[j] == -1) {
    match_vector[j] <- "Human↑ / Mouse↓"
  } else if (h_dir[j] == -1 & m_dir[j] == 1) {
    match_vector[j] <- "Human↓ / Mouse↑"
  } else if (h_dir[j] == 0 & m_dir[j] != 0) {
    match_vector[j] <- "Human↔ / Mouse≠"
  } else if (h_dir[j] != 0 & m_dir[j] == 0) {
    match_vector[j] <- "Human≠ / Mouse↔"
  } else {
    match_vector[j] <- "No Change"
  }
}

  # Save overall match score
direction_match$match_score[i] <- sum(match_vector == "Match")

  # Save transition-level matches
  direction_match[i, transition_labels] <- match_vector
  

 # Determine regulation trend based on 0 vs 96 timepoint
direction_match$human_trend[i] <- if (h_expr[length(h_expr)] > h_expr[1]) {
  "Upregulated"
} else if (h_expr[length(h_expr)] < h_expr[1]) {
  "Downregulated"
} else {
  "No Change"
}

direction_match$mouse_trend[i] <- if (m_expr[length(m_expr)] > m_expr[1]) {
  "Upregulated"
} else if (m_expr[length(m_expr)] < m_expr[1]) {
  "Downregulated"
} else {
  "No Change"
}
}

# Add percentage of matching directions
direction_match$percent_match <- round(direction_match$match_score / direction_match$total_transitions * 100, 2)


# Final safeguard to remove any unintended columns (e.g., "T0.0_to_T12.1")
expected_cols <- c(
  "gene_name", "match_score", "total_transitions", transition_labels,
  "human_trend", "mouse_trend", "percent_match"
)
direction_match <- direction_match[, intersect(expected_cols, colnames(direction_match))]

# View result
print(direction_match)
```

    ##    gene_name match_score total_transitions       T0_to_T12      T12_to_T24
    ## 1      abcc2           1                 4 Human↑ / Mouse↓           Match
    ## 2      ptprn           2                 4           Match Human↑ / Mouse↓
    ## 3      fgfr3           3                 4           Match Human↑ / Mouse↓
    ## 4       lhx5           2                 4           Match Human↑ / Mouse↓
    ## 5        bik           1                 4 Human↑ / Mouse↓           Match
    ## 6        gck           4                 4           Match           Match
    ## 7       oas3           1                 4 Human↑ / Mouse↓ Human↑ / Mouse↓
    ## 8      rbm24           2                 4 Human↑ / Mouse↓           Match
    ## 9       bche           3                 4           Match           Match
    ## 10   ppp1r3c           1                 4 Human↑ / Mouse↓ Human↑ / Mouse↓
    ## 11      psd4           3                 4           Match Human↑ / Mouse↓
    ## 12     chac1           2                 4 Human↓ / Mouse↑           Match
    ## 13    atp1a4           3                 4           Match           Match
    ## 14       cpm           3                 4 Human↑ / Mouse↓           Match
    ## 15      ier3           2                 4           Match Human↑ / Mouse↓
    ## 16   slc40a1           2                 4 Human↑ / Mouse↓ Human↓ / Mouse↑
    ## 17    cyp1a1           4                 4           Match           Match
    ## 18      rhcg           2                 4           Match Human↓ / Mouse↑
    ## 19     padi3           0                 4 Human↑ / Mouse↓ Human↑ / Mouse↓
    ## 20    itga10           2                 4           Match           Match
    ## 21     gadl1           0                 4 Human↓ / Mouse↑ Human↓ / Mouse↑
    ## 22     ccdc3           1                 4 Human↑ / Mouse↓ Human↑ / Mouse↓
    ## 23   pstpip2           0                 4 Human↑ / Mouse↓ Human↑ / Mouse↓
    ## 24     robo4           2                 4 Human↓ / Mouse↑           Match
    ## 25     htr5a           3                 4           Match           Match
    ## 26      cbr3           2                 4           Match Human↓ / Mouse↑
    ## 27    elavl4           1                 4 Human↑ / Mouse↓           Match
    ## 28      otx2           2                 4           Match Human↑ / Mouse↓
    ## 29      mei1           2                 4 Human↑ / Mouse↓ Human↑ / Mouse↓
    ## 30      igf2           1                 4 Human↑ / Mouse↓           Match
    ## 31     gap43           2                 4 Human↓ / Mouse↑           Match
    ## 32      etv4           3                 4 Human↓ / Mouse↑           Match
    ## 33      ly6h           1                 4 Human↑ / Mouse↓ Human↑ / Mouse↓
    ## 34     dppa3           1                 4           Match Human↑ / Mouse↓
    ## 35      dll1           1                 4           Match Human↑ / Mouse↓
    ## 36       sp5           2                 4 Human↑ / Mouse↓           Match
    ##         T24_to_T48      T48_to_T96   human_trend   mouse_trend percent_match
    ## 1  Human↑ / Mouse↓ Human↑ / Mouse↓ Downregulated Downregulated            25
    ## 2  Human↓ / Mouse↑           Match Downregulated Downregulated            50
    ## 3            Match           Match   Upregulated Downregulated            75
    ## 4            Match Human↓ / Mouse↑   Upregulated Downregulated            50
    ## 5  Human↑ / Mouse↓ Human↓ / Mouse↑   Upregulated Downregulated            25
    ## 6            Match           Match Downregulated Downregulated           100
    ## 7            Match Human↓ / Mouse↑   Upregulated Downregulated            25
    ## 8  Human↑ / Mouse↓           Match   Upregulated Downregulated            50
    ## 9            Match Human↑ / Mouse↓   Upregulated   Upregulated            75
    ## 10 Human↓ / Mouse↑           Match   Upregulated Downregulated            25
    ## 11           Match           Match   Upregulated Downregulated            75
    ## 12           Match Human↓ / Mouse↑ Downregulated Downregulated            50
    ## 13 Human↑ / Mouse↓           Match   Upregulated Downregulated            75
    ## 14           Match           Match Downregulated Downregulated            75
    ## 15 Human↓ / Mouse↑           Match Downregulated Downregulated            50
    ## 16           Match           Match   Upregulated   Upregulated            50
    ## 17           Match           Match   Upregulated   Upregulated           100
    ## 18 Human↓ / Mouse↑           Match   Upregulated   Upregulated            50
    ## 19 Human↓ / Mouse↑ Human↓ / Mouse↑   Upregulated Downregulated             0
    ## 20 Human↑ / Mouse↓ Human↓ / Mouse↑   Upregulated Downregulated            50
    ## 21 Human↓ / Mouse↑ Human↑ / Mouse↓ Downregulated Downregulated             0
    ## 22           Match Human↓ / Mouse↑   Upregulated Downregulated            25
    ## 23 Human↓ / Mouse↑ Human↑ / Mouse↓   Upregulated Downregulated             0
    ## 24 Human↑ / Mouse↓           Match   Upregulated   Upregulated            50
    ## 25           Match Human↑ / Mouse↓ Downregulated Downregulated            75
    ## 26 Human↓ / Mouse↑           Match Downregulated   Upregulated            50
    ## 27 Human↑ / Mouse↓ Human↓ / Mouse↑   Upregulated   Upregulated            25
    ## 28           Match Human↓ / Mouse↑   Upregulated Downregulated            50
    ## 29           Match           Match   Upregulated Downregulated            50
    ## 30 Human↓ / Mouse↑ Human↑ / Mouse↓   Upregulated   Upregulated            25
    ## 31           Match Human↓ / Mouse↑ Downregulated Downregulated            50
    ## 32           Match           Match Downregulated Downregulated            75
    ## 33           Match Human↓ / Mouse↑   Upregulated Downregulated            25
    ## 34 Human↓ / Mouse↑ Human↑ / Mouse↓ Downregulated Downregulated            25
    ## 35 Human↑ / Mouse↓ Human↓ / Mouse↑   Upregulated Downregulated            25
    ## 36           Match Human↑ / Mouse↓   Upregulated   Upregulated            50

The most notable data we see from this table is that Gck and Cyp1a1
expression transcripts in human and mouse matched at all 4 time points
indicating that when the human gene is upregulated/downregulated, the
mouse version of the gene saw the same trend. When comparing this to the
line plots, GCK and Cyp1a1 showed consistent trends in both the
directional match data and the line plots. Another notable difference is
the matching score of 50% for IER3. This is surprising because when
comparing the human and mouse plots they appear to be very similar and
you would expect a higher matching score.

Now we will run similar code to above but this time we will calculate
the log fold change of TPM to see if there will be any differences with
the data

``` r
# Initialize new results table
lfc_match <- data.frame(
  gene_name = common_genes,
  match_score = NA_integer_,
  total_transitions = n_transitions,
  matrix(NA_character_, nrow = length(common_genes), ncol = n_transitions),
  human_trend = NA_character_,
  mouse_trend = NA_character_,
  percent_match = NA_real_
)

colnames(lfc_match)[4:(3 + n_transitions)] <- transition_labels

# Loop through each gene
for (i in seq_along(common_genes)) {
  gene <- common_genes[i]

  h_row <- filtered_human_tpm[tolower(filtered_human_tpm$gene_name) == gene, ]
  m_row <- mouse_conserved_tpm[tolower(mouse_conserved_tpm$gene_name) == gene, ]

  h_expr <- sapply(time_points, function(tp) {
    rowMeans(h_row[, grep(paste0("gfp_", tp, "_"), colnames(h_row)), drop = FALSE])
  })
  m_expr <- sapply(time_points, function(tp) {
    rowMeans(m_row[, grep(paste0("WT_", tp, "_"), colnames(m_row)), drop = FALSE])
  })

  # Add pseudocount to avoid log(0)
  h_lfc <- log2(h_expr[-1] + 1) - log2(h_expr[-length(h_expr)] + 1)
  m_lfc <- log2(m_expr[-1] + 1) - log2(m_expr[-length(m_expr)] + 1)

  # Classify fold change direction
  match_vector <- character(length = n_transitions)
  for (j in seq_len(n_transitions)) {
    if (sign(h_lfc[j]) == sign(m_lfc[j])) {
      match_vector[j] <- "Match"
    } else if (h_lfc[j] > 0 & m_lfc[j] < 0) {
      match_vector[j] <- "Human↑ / Mouse↓"
    } else if (h_lfc[j] < 0 & m_lfc[j] > 0) {
      match_vector[j] <- "Human↓ / Mouse↑"
    } else if (h_lfc[j] == 0 & m_lfc[j] != 0) {
      match_vector[j] <- "Human↔ / Mouse≠"
    } else if (h_lfc[j] != 0 & m_lfc[j] == 0) {
      match_vector[j] <- "Human≠ / Mouse↔"
    } else {
      match_vector[j] <- "No Change"
    }
  }

  # Store results
  lfc_match$match_score[i] <- sum(match_vector == "Match")
  lfc_match[i, transition_labels] <- match_vector

  # Overall trend based on 0h to 96h LFC
  overall_h_lfc <- log2(h_expr[length(h_expr)] + 1) - log2(h_expr[1] + 1)
  overall_m_lfc <- log2(m_expr[length(m_expr)] + 1) - log2(m_expr[1] + 1)

  lfc_match$human_trend[i] <- if (overall_h_lfc > 0) {
    "Upregulated"
  } else if (overall_h_lfc < 0) {
    "Downregulated"
  } else {
    "No Change"
  }

  lfc_match$mouse_trend[i] <- if (overall_m_lfc > 0) {
    "Upregulated"
  } else if (overall_m_lfc < 0) {
    "Downregulated"
  } else {
    "No Change"
  }

  lfc_match$percent_match[i] <- round(lfc_match$match_score[i] / n_transitions * 100, 2)
}

#Optional if you would like to save the directional matching tables in an easier format to read, use the code below to save the tables in an excel sheet. Make sure to install openxlsx and load in the library. Otherwise the code below will not work. 

# Saving the direction_match and lfc_match as an excel spreadsheet 

# Output directory and file
output_dir <- "/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/brainiacs/brainiacs_dox_genes/Tables"
output_file <- file.path(output_dir, "gene_expression_matches_pretty.xlsx")

# Create workbook
wb <- createWorkbook()

# Formatting style
header_style <- createStyle(
  fontSize = 12, fontColour = "#FFFFFF", fgFill = "#4F81BD",
  halign = "center", textDecoration = "bold", border = "Bottom"
)

# Write direction_match
addWorksheet(wb, "Direction_Match")
writeData(wb, "Direction_Match", direction_match, headerStyle = header_style)
addStyle(wb, "Direction_Match", header_style, rows = 1, cols = 1:ncol(direction_match), gridExpand = TRUE)
freezePane(wb, "Direction_Match", firstRow = TRUE)
setColWidths(wb, "Direction_Match", cols = 1:ncol(direction_match), widths = "auto")
addFilter(wb, "Direction_Match", rows = 1, cols = 1:ncol(direction_match))

# Write lfc_match
addWorksheet(wb, "LFC_Match")
writeData(wb, "LFC_Match", lfc_match, headerStyle = header_style)
addStyle(wb, "LFC_Match", header_style, rows = 1, cols = 1:ncol(lfc_match), gridExpand = TRUE)
freezePane(wb, "LFC_Match", firstRow = TRUE)
setColWidths(wb, "LFC_Match", cols = 1:ncol(lfc_match), widths = "auto")
addFilter(wb, "LFC_Match", rows = 1, cols = 1:ncol(lfc_match))

# Save
saveWorkbook(wb, output_file, overwrite = TRUE)
```

We see that there is no difference between the direction (TPM) and log
fold change tables and that the trends seen in one, are seen in the
other.

Next we will construct a bar plot to see how many of the genes in Mouse
and Human are primarily being downregulated or upregulated across the
time points.

``` r
# Make sure lfc_match is a proper data frame
lfc_match <- as.data.frame(lfc_match)

# Pivot longer to combine human and mouse trends
trend_counts <- lfc_match %>%
  dplyr::select(gene_name, human_trend, mouse_trend) %>%
  pivot_longer(
    cols = c(human_trend, mouse_trend),
    names_to = "Species",
    values_to = "Trend"
  ) %>%
  mutate(
    Species = recode(Species,
                     human_trend = "Human",
                     mouse_trend = "Mouse")
  ) %>%
  group_by(Species, Trend) %>%
  summarise(n = n(), .groups = "drop")

# Plot the trend counts
ggplot(trend_counts, aes(x = Trend, y = n, fill = Species)) +
  geom_col(position = "dodge") +
  labs(
    title = "Gene Regulation Trends by Species (LFC 0h -> 96h)",
    x = "Regulation Trend",
    y = "Number of Genes"
  ) +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  theme(text = element_text(size = 12))
```

![](Final_files/figure-gfm/Plotting%20Bar%20Plot-1.png)<!-- --> We see
from the bar plot that from the 0 time point to the 96 hour time, there
is more genes downregulated in Mouse and more genes upregulated in Human
after being treated with Doxycycline.

## KEGG Pathway Annotation of Conserved Orthologues

We found 36 total significant genes that are common between human and
mouse after doxycyclene treatment. We utilized the KEGG database to
build a table of the genes, their Entrez ID, and the descriptions of the
pathways the gene is affiliated with.

``` r
library(biomaRt)

# Set up mouse and human marts
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# mouse gene list
mouse_list <- c("Htr5a", "Ier3", "Cbr3", "Ccdc3", "Dppa3", "Mei1", "Rbm24", "Pstpip2", "Dll1", "Bche", "Cpm", "Otx2", "Slc40a1", "Bik", "Ptprn", "Gck", "Ppp1r3c", "Elavl4", "Oas3", "Abcc2", "Lhx5", "Padi3", "Psd4", "Gadl1", "Fgfr3", "Cyp1a1", "Chac1", "Etv4","Rhcg", "Sp5", "Robo4", "Itga10", "Ly6h", "Gap43", "Igf2", "Atp1a4")

#Builds entrez gene list
mouse_ids <- getBM(attributes = c("mgi_symbol", "entrezgene_id"),
                   filters = "mgi_symbol",
                   values = mouse_list,
                   mart = mouse)

# Now for human
human_ids <- data.frame(
  hgnc_symbol = c("Htr5a", "Ier3", "Cbr3", "Ccdc3", "Dppa3", "Mei1", "Rbm24", "Pstpip2",
                  "Dll1", "Bche", "Cpm", "Otx2", "Slc40a1", "Bik", "Ptprn", "Gck", "Ppp1r3c",
                  "Elavl4", "Oas3", "Abcc2", "Lhx5", "Padi3", "Psd4", "Gadl1", "Fgfr3",
                  "Cyp1a1", "Chac1", "Etv4", "Rhcg", "Sp5", "Robo4", "Itga10", "Ly6h",
                  "Gap43", "Igf2", "Atp1a4"),
  human_entrez_id = c(3361, 8870, 874, 83643, 359787, 150365, 221662, 9050,
                28514, 590, 1368, 5015, 30061, 638, 5798, 2645, 5509,
                1995, 4940, 1244, 5589, 51702, 26037, 144165, 2261,
                1543, 79094, 2118, 51499, 222489, 54538, 8515, 4061,
                2596, 3481, 479)
)
merged_ids <- merge(mouse_ids, human_ids, 
                    by.x = "mgi_symbol", by.y = "hgnc_symbol", 
                    all.x = TRUE)


# Rename for clarity
colnames(merged_ids) <- c("Gene", "Mouse_Entrez_ID", "Human_Entrez_ID")

# Function to get KEGG pathways
get_kegg_pathways <- function(entrez_id, species_prefix = "mmu") {
  if (is.na(entrez_id) || entrez_id == "") return(NA)
  url <- paste0("https://rest.kegg.jp/get/", species_prefix, ":", entrez_id)
  response <- tryCatch(readLines(url, warn = FALSE), error = function(e) return(character(0)))
  if (length(response) == 0) return(NA)
  
  pathways <- grep("^PATHWAY", response, value = TRUE)
  if (length(pathways) == 0) return(NA)
  gsub("PATHWAY\\s+", "", pathways)
}

# Apply KEGG retrieval to each row
merged_ids$Mouse_KEGG_Pathways <- sapply(merged_ids$Mouse_Entrez_ID, function(id) {
  pw <- get_kegg_pathways(id, species_prefix = "mmu")
  if (all(is.na(pw))) return(NA)
  paste(pw, collapse = "; ")
})

merged_ids$Human_KEGG_Pathways <- sapply(merged_ids$Human_Entrez_ID, function(id) {
  pw <- get_kegg_pathways(id, species_prefix = "hsa")
  if (all(is.na(pw))) return(NA)
  paste(pw, collapse = "; ")
})

# Final table
kegg_pathway_table <- merged_ids
print(kegg_pathway_table[, c("Gene", "Mouse_KEGG_Pathways", "Human_KEGG_Pathways")])
```

    ##       Gene                                                Mouse_KEGG_Pathways
    ## 1    Abcc2                                    mmu01523  Antifolate resistance
    ## 2   Atp1a4                               mmu04022  cGMP-PKG signaling pathway
    ## 3     Bche                                                               <NA>
    ## 4      Bik                                     mmu01522  Endocrine resistance
    ## 5     Cbr3                              mmu00590  Arachidonic acid metabolism
    ## 6    Ccdc3                                                               <NA>
    ## 7    Chac1                                   mmu00480  Glutathione metabolism
    ## 8      Cpm                                                               <NA>
    ## 9   Cyp1a1                             mmu00140  Steroid hormone biosynthesis
    ## 10    Dll1                                     mmu01522  Endocrine resistance
    ## 11   Dppa3                                                               <NA>
    ## 12  Elavl4                                                               <NA>
    ## 13    Etv4                  mmu05202  Transcriptional misregulation in cancer
    ## 14   Fgfr3                mmu01521  EGFR tyrosine kinase inhibitor resistance
    ## 15   Gadl1                                  mmu00410  beta-Alanine metabolism
    ## 16   Gap43                                                               <NA>
    ## 17     Gck                             mmu00010  Glycolysis / Gluconeogenesis
    ## 18   Htr5a                                mmu04020  Calcium signaling pathway
    ## 19    Ier3                                                               <NA>
    ## 20    Igf2                                   mmu04010  MAPK signaling pathway
    ## 21  Itga10                               mmu04151  PI3K-Akt signaling pathway
    ## 22    Lhx5 mmu04550  Signaling pathways regulating pluripotency of stem cells
    ## 23    Ly6h                                                               <NA>
    ## 24    Mei1                                                               <NA>
    ## 25    Oas3                      mmu04621  NOD-like receptor signaling pathway
    ## 26    Otx2                                                               <NA>
    ## 27   Padi3                                                               <NA>
    ## 28 Ppp1r3c                                mmu04910  Insulin signaling pathway
    ## 29    Psd4                                              mmu04144  Endocytosis
    ## 30 Pstpip2                                                               <NA>
    ## 31   Ptprn                                 mmu04940  Type I diabetes mellitus
    ## 32   Rbm24                                                               <NA>
    ## 33    Rhcg                                                               <NA>
    ## 34   Robo4                                                               <NA>
    ## 35 Slc40a1                                        mmu04081  Hormone signaling
    ## 36     Sp5                                                               <NA>
    ##                                      Human_KEGG_Pathways
    ## 1                        hsa01523  Antifolate resistance
    ## 2                    hsa00190  Oxidative phosphorylation
    ## 3                                                   <NA>
    ## 4                         hsa01522  Endocrine resistance
    ## 5                  hsa00590  Arachidonic acid metabolism
    ## 6                                                   <NA>
    ## 7                       hsa00480  Glutathione metabolism
    ## 8                                                   <NA>
    ## 9                 hsa00140  Steroid hormone biosynthesis
    ## 10                        hsa01522  Endocrine resistance
    ## 11                                                  <NA>
    ## 12                                                  <NA>
    ## 13     hsa05202  Transcriptional misregulation in cancer
    ## 14   hsa01521  EGFR tyrosine kinase inhibitor resistance
    ## 15                       hsa04310  Wnt signaling pathway
    ## 16                                                  <NA>
    ## 17                hsa00010  Glycolysis / Gluconeogenesis
    ## 18                   hsa04020  Calcium signaling pathway
    ## 19                                                  <NA>
    ## 20                      hsa04010  MAPK signaling pathway
    ## 21                  hsa04151  PI3K-Akt signaling pathway
    ## 22 hsa04141  Protein processing in endoplasmic reticulum
    ## 23                                                  <NA>
    ## 24                                                  <NA>
    ## 25         hsa04621  NOD-like receptor signaling pathway
    ## 26                                                  <NA>
    ## 27                                                  <NA>
    ## 28                   hsa04910  Insulin signaling pathway
    ## 29                      hsa04015  Rap1 signaling pathway
    ## 30                                                  <NA>
    ## 31                    hsa04940  Type I diabetes mellitus
    ## 32                                                  <NA>
    ## 33                                                  <NA>
    ## 34                                                  <NA>
    ## 35                           hsa04081  Hormone signaling
    ## 36                                                  <NA>

The above image shows a table of the 36 human and mouse orthologous
genes, their Entrez IDs, and the KEGG pathway annotations. From the KEGG
database, we found that many genes between human and mouse had similar
pathways. Common pathways highlighted key biological processes,
including signal transduction, metabolic pathways, and
disease-associated mechanisms. NA values indicate genes without an
annotated KEGG pathway in the respective species.

## Gene Expression Changes Over Time (RNA-seq IGV Snapshots)

For each gene, human RNA-seq IGV visualizations are shown on the left,
and mouse RNA-seq IGV visualizations are shown on the right.

RNA-seq IGV Visualization of HTR5A Across 0–96 Hour Time Points \| HTR5A
\| \| \|———————–\|———————–\| \|
![](RNASEQ_Mouse_and_Human_IGV/igv_snapshot-HTR5A-1-Human.png) \|
![](RNASEQ_Mouse_and_Human_IGV/igv_snapshot-Htr.png) \| \|
![](RNASEQ_Mouse_and_Human_IGV/igv_snapshot-HTR5A-2-Human.png) \|
![](RNASEQ_Mouse_and_Human_IGV/igv_snapshot-Htr-2.png) \| \|
![](RNASEQ_Mouse_and_Human_IGV/igv_snapshot-HTR5A-3-Human.png) \|
![](RNASEQ_Mouse_and_Human_IGV/igv_snapshot-Htr-3.png) \|

RNA-seq IGV Visualization of IER3 Across 0–96 Hour Time Points \| IER3
\| \| \|———————–\|———————–\| \|
![](RNASEQ_Mouse_and_Human_IGV/igv_snapshot-IER3-1-Human.png) \|
![](RNASEQ_Mouse_and_Human_IGV/igv_snapshot-Ier3.png) \| \|
![](RNASEQ_Mouse_and_Human_IGV/igv_snapshot-IER3-2-Human.png) \|
![](RNASEQ_Mouse_and_Human_IGV/igv_snapshot-ler3-2.png) \| \|
![](RNASEQ_Mouse_and_Human_IGV/igv_snapshot-IER3-3-Human.png) \|
![](RNASEQ_Mouse_and_Human_IGV/igv_snapshot-ler3-3.png) \|

## Chromatin Accessibility Profiles Across Time (ATAC-seq IGV Snapshots)

For each gene, human ATAC-seq IGV visualizations are shown on the left,
and mouse ATAC-seq IGV visualizations are shown on the right. Human IGV
images display chromatin accessibility patterns across 0 to 2.5 hours
following doxycycline treatment, with three biological replicates per
time point (0h, 0.5h, 1h, 1.5h, 2h, and 2.5h). In contrast, mouse IGV
images display chromatin accessibility across 0 to 2.5 hours by
comparing wild-type (WT) and knockout (KO) conditions at 0, 30, 60, 90,
120, and 150 minutes.

| ABCC2                                                     |                                                           |
|-----------------------------------------------------------|-----------------------------------------------------------|
| ![](ATAC_SEQ_Mouse_and_Human_IGV/human_atacseq_ABCC2.png) | ![](ATAC_SEQ_Mouse_and_Human_IGV/mouse_atacseq_ABCC2.png) |

ATAC-seq IGV visualization of ABCC2 in human and mouse

| ATP1A4                                                     |                                                            |
|------------------------------------------------------------|------------------------------------------------------------|
| ![](ATAC_SEQ_Mouse_and_Human_IGV/human_atacseq_ATP1A4.png) | ![](ATAC_SEQ_Mouse_and_Human_IGV/mouse_atacseq_ATP1A4.png) |

ATAC-seq IGV visualization of ATP1A4 in human and mouse

| BCHE                                                     |                                                          |
|----------------------------------------------------------|----------------------------------------------------------|
| ![](ATAC_SEQ_Mouse_and_Human_IGV/human_atacseq_BCHE.png) | ![](ATAC_SEQ_Mouse_and_Human_IGV/mouse_atacseq_BCHE.png) |

ATAC-seq IGV visualization of BCHE in human and mouse

| BIK                                                     |                                                         |
|---------------------------------------------------------|---------------------------------------------------------|
| ![](ATAC_SEQ_Mouse_and_Human_IGV/human_atacseq_BIK.png) | ![](ATAC_SEQ_Mouse_and_Human_IGV/mouse_atacseq_BIK.png) |

ATAC-seq IGV visualization of BIK in human and mouse

| CBR3                                                     |                                                          |
|----------------------------------------------------------|----------------------------------------------------------|
| ![](ATAC_SEQ_Mouse_and_Human_IGV/human_atacseq_CBR3.png) | ![](ATAC_SEQ_Mouse_and_Human_IGV/mouse_atacseq_CBR3.png) |

ATAC-seq IGV visualization of CBR3 in human and mouse

| CCDC3                                                     |                                                           |
|-----------------------------------------------------------|-----------------------------------------------------------|
| ![](ATAC_SEQ_Mouse_and_Human_IGV/human_atacseq_CCDC3.png) | ![](ATAC_SEQ_Mouse_and_Human_IGV/mouse_atacseq_CCDC3.png) |

ATAC-seq IGV visualization of CCDC3 in human and mouse

| CHAC1                                                     |                                                           |
|-----------------------------------------------------------|-----------------------------------------------------------|
| ![](ATAC_SEQ_Mouse_and_Human_IGV/human_atacseq_CHAC1.png) | ![](ATAC_SEQ_Mouse_and_Human_IGV/mouse_atacseq_CHAC1.png) |

ATAC-seq IGV visualization of CHAC1 in human and mouse

| CPM                                                     |                                                         |
|---------------------------------------------------------|---------------------------------------------------------|
| ![](ATAC_SEQ_Mouse_and_Human_IGV/human_atacseq_CPM.png) | ![](ATAC_SEQ_Mouse_and_Human_IGV/mouse_atacseq_CPM.png) |

ATAC-seq IGV visualization of CPM in human and mouse

| DLL1                                                     |                                                          |
|----------------------------------------------------------|----------------------------------------------------------|
| ![](ATAC_SEQ_Mouse_and_Human_IGV/human_atacseq_DLL1.png) | ![](ATAC_SEQ_Mouse_and_Human_IGV/mouse_atacseq_DLL1.png) |

ATAC-seq IGV visualization of DLL1 in human and mouse

| DPPA3                                                     |                                                           |
|-----------------------------------------------------------|-----------------------------------------------------------|
| ![](ATAC_SEQ_Mouse_and_Human_IGV/human_atacseq_DPPA3.png) | ![](ATAC_SEQ_Mouse_and_Human_IGV/mouse_atacseq_DPPA3.png) |

ATAC-seq IGV visualization of DPPA3 in human and mouse

| ELAVL4                                                     |                                                            |
|------------------------------------------------------------|------------------------------------------------------------|
| ![](ATAC_SEQ_Mouse_and_Human_IGV/human_atacseq_ELAVL4.png) | ![](ATAC_SEQ_Mouse_and_Human_IGV/mouse_atacseq_ELAVL4.png) |

ATAC-seq IGV visualization of ELAVL4 in human and mouse

| ETV4                                                     |                                                          |
|----------------------------------------------------------|----------------------------------------------------------|
| ![](ATAC_SEQ_Mouse_and_Human_IGV/human_atacseq_ETV4.png) | ![](ATAC_SEQ_Mouse_and_Human_IGV/mouse_atacseq_ETV4.png) |

ATAC-seq IGV visualization of ETV4 in human and mouse

| FGFR3                                                     |                                                           |
|-----------------------------------------------------------|-----------------------------------------------------------|
| ![](ATAC_SEQ_Mouse_and_Human_IGV/human_atacseq_FGFR3.png) | ![](ATAC_SEQ_Mouse_and_Human_IGV/mouse_atacseq_FGFR3.png) |

ATAC-seq IGV visualization of FGFR3 in human and mouse

| GADL1                                                     |                                                           |
|-----------------------------------------------------------|-----------------------------------------------------------|
| ![](ATAC_SEQ_Mouse_and_Human_IGV/human_atacseq_GADL1.png) | ![](ATAC_SEQ_Mouse_and_Human_IGV/mouse_atacseq_GADL1.png) |

ATAC-seq IGV visualization of GADL1 in human and mouse

| GAP43                                                     |                                                           |
|-----------------------------------------------------------|-----------------------------------------------------------|
| ![](ATAC_SEQ_Mouse_and_Human_IGV/human_atacseq_GAP43.png) | ![](ATAC_SEQ_Mouse_and_Human_IGV/mouse_atacseq_GAP43.png) |

ATAC-seq IGV visualization of GAP43 in human and mouse

| GCK                                                     |                                                         |
|---------------------------------------------------------|---------------------------------------------------------|
| ![](ATAC_SEQ_Mouse_and_Human_IGV/human_atacseq_GCK.png) | ![](ATAC_SEQ_Mouse_and_Human_IGV/mouse_atacseq_GCK.png) |

ATAC-seq IGV visualization of GCK in human and mouse

| HTR5A                                                     |                                                           |
|-----------------------------------------------------------|-----------------------------------------------------------|
| ![](ATAC_SEQ_Mouse_and_Human_IGV/human_atacseq_HTR5A.png) | ![](ATAC_SEQ_Mouse_and_Human_IGV/mouse_atacseq_HTR5A.png) |

ATAC-seq IGV visualization of HTR5A in human and mouse

| IER3                                                     |                                                          |
|----------------------------------------------------------|----------------------------------------------------------|
| ![](ATAC_SEQ_Mouse_and_Human_IGV/human_atacseq_IER3.png) | ![](ATAC_SEQ_Mouse_and_Human_IGV/mouse_atacseq_IER3.png) |

ATAC-seq IGV visualization of IER3 in human and mouse

| IGF2                                                     |                                                          |
|----------------------------------------------------------|----------------------------------------------------------|
| ![](ATAC_SEQ_Mouse_and_Human_IGV/human_atacseq_IGF2.png) | ![](ATAC_SEQ_Mouse_and_Human_IGV/mouse_atacseq_IGF2.png) |

ATAC-seq IGV visualization of IGF2 in human and mouse

| ITGA10                                                     |                                                            |
|------------------------------------------------------------|------------------------------------------------------------|
| ![](ATAC_SEQ_Mouse_and_Human_IGV/human_atacseq_ITGA10.png) | ![](ATAC_SEQ_Mouse_and_Human_IGV/mouse_atacseq_ITGA10.png) |

ATAC-seq IGV visualization of ITGA10 in human and mouse

| LHX5                                                     |                                                          |
|----------------------------------------------------------|----------------------------------------------------------|
| ![](ATAC_SEQ_Mouse_and_Human_IGV/human_atacseq_LHX5.png) | ![](ATAC_SEQ_Mouse_and_Human_IGV/mouse_atacseq_LHX5.png) |

ATAC-seq IGV visualization of LHX5 in human and mouse

| LY6H                                                     |                                                          |
|----------------------------------------------------------|----------------------------------------------------------|
| ![](ATAC_SEQ_Mouse_and_Human_IGV/human_atacseq_LY6H.png) | ![](ATAC_SEQ_Mouse_and_Human_IGV/mouse_atacseq_LY6H.png) |

ATAC-seq IGV visualization of LY6H in human and mouse

| MEI1                                                     |                                                          |
|----------------------------------------------------------|----------------------------------------------------------|
| ![](ATAC_SEQ_Mouse_and_Human_IGV/human_atacseq_MEI1.png) | ![](ATAC_SEQ_Mouse_and_Human_IGV/mouse_atacseq_MEI1.png) |

ATAC-seq IGV visualization of MEI1 in human and mouse

| OAS3                                                     |                                                          |
|----------------------------------------------------------|----------------------------------------------------------|
| ![](ATAC_SEQ_Mouse_and_Human_IGV/human_atacseq_OAS3.png) | ![](ATAC_SEQ_Mouse_and_Human_IGV/mouse_atacseq_OAS3.png) |

ATAC-seq IGV visualization of OAS3 in human and mouse

| OTX2                                                     |                                                          |
|----------------------------------------------------------|----------------------------------------------------------|
| ![](ATAC_SEQ_Mouse_and_Human_IGV/human_atacseq_OTX2.png) | ![](ATAC_SEQ_Mouse_and_Human_IGV/mouse_atacseq_OTX2.png) |

ATAC-seq IGV visualization of OTX2 in human and mouse

| PADI3                                                     |                                                           |
|-----------------------------------------------------------|-----------------------------------------------------------|
| ![](ATAC_SEQ_Mouse_and_Human_IGV/human_atacseq_PADI3.png) | ![](ATAC_SEQ_Mouse_and_Human_IGV/mouse_atacseq_PADI3.png) |

ATAC-seq IGV visualization of PADI3 in human and mouse

| PPP1R3C                                                     |                                                             |
|-------------------------------------------------------------|-------------------------------------------------------------|
| ![](ATAC_SEQ_Mouse_and_Human_IGV/human_atacseq_PPP1R3C.png) | ![](ATAC_SEQ_Mouse_and_Human_IGV/mouse_atacseq_PPP1R3C.png) |

ATAC-seq IGV visualization of PPP1R3C in human and mouse

| RBM24                                                     |                                                           |
|-----------------------------------------------------------|-----------------------------------------------------------|
| ![](ATAC_SEQ_Mouse_and_Human_IGV/human_atacseq_RBM24.png) | ![](ATAC_SEQ_Mouse_and_Human_IGV/mouse_atacseq_RBM24.png) |

ATAC-seq IGV visualization of RBM24 in human and mouse

| RHCG                                                     |                                                          |
|----------------------------------------------------------|----------------------------------------------------------|
| ![](ATAC_SEQ_Mouse_and_Human_IGV/human_atacseq_RHCG.png) | ![](ATAC_SEQ_Mouse_and_Human_IGV/mouse_atacseq_RHCG.png) |

ATAC-seq IGV visualization of RHCG in human and mouse

| ROBO4                                                     |                                                           |
|-----------------------------------------------------------|-----------------------------------------------------------|
| ![](ATAC_SEQ_Mouse_and_Human_IGV/human_atacseq_ROBO4.png) | ![](ATAC_SEQ_Mouse_and_Human_IGV/mouse_atacseq_ROBO4.png) |

ATAC-seq IGV visualization of ROBO4 in human and mouse

| SLC40A1                                                     |                                                             |
|-------------------------------------------------------------|-------------------------------------------------------------|
| ![](ATAC_SEQ_Mouse_and_Human_IGV/human_atacseq_SLC40A1.png) | ![](ATAC_SEQ_Mouse_and_Human_IGV/mouse_atacseq_SLC40A1.png) |

ATAC-seq IGV visualization of SLC40A1 in human and mouse

| SP5                                                     |                                                         |
|---------------------------------------------------------|---------------------------------------------------------|
| ![](ATAC_SEQ_Mouse_and_Human_IGV/human_atacseq_SP5.png) | ![](ATAC_SEQ_Mouse_and_Human_IGV/mouse_atacseq_SP5.png) |

ATAC-seq IGV visualization of SP5 in human and mouse

The IGV plots of the ATAC-seq peak signals revealed minimal changes
across time points for all 36 genes in the human dataset, with genes
DLL1, IER3, LY6H, and PADI3 showing no detectable peaks at any time
point. In contrast, the mouse ATAC-seq data showed a more noticeable
pattern, with most genes exhibiting their highest peaks at the 90-minute
knockout (KO) time point, suggesting increased chromatin accessibility
specifically in response to the KO. Notably, the gene PSTPIP2 was the
only one in the mouse dataset that showed no peaks across all time
points. Interestingly, for each of the genes that lacked peaks in one
species, there were clear peaks in the other, PSTPIP2 had peaks in
human, while DLL1, IER3, LY6H, and PADI3 had peaks in mouse,
highlighting a divergence in chromatin accessibility between species for
specific genes, which may contribute to differences in gene regulation.

## Conclusion

Our analysis identified 36 orthologous genes exhibiting significant
expression changes across time in both human and mouse models after
doxycycline treatment. Analysis of regulation showed that about 50% of
genes had conserved up or down regulation. Log fold change analysis
showed that changes were consistent between raw data and
transformations. RNA_seq IGV plots supported observed gene expression
trends, while ATAC_seq IGV plots showed increased chromatin
accessibility. Overall, we conducted a comprehensive comparative
analysis of gene expression dynamics, and human and mouse
transcriptional responses to doxycycline were broadly conserved.
