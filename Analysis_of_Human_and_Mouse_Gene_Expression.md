Analysis_of_Human_and_Mouse_Gene_Expressions
================

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

![](Analysis_of_Human_and_Mouse_Gene_Expression_files/figure-gfm/Plotting%20TPM%20Data%20from%20Mouse%20and%20Human%20Data-1.png)<!-- -->

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

![](Analysis_of_Human_and_Mouse_Gene_Expression_files/figure-gfm/Plotting%20TPM%20Data%20from%20Mouse%20and%20Human%20Data-2.png)<!-- -->

\# CREATING DIRECTION VECTORS FROM THE HUMAN AND MOUSE DATA TO SEE IF
THEIR TRANSCRIPTS ARE GOING IN THE SAME DIRECTIONS

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

# Now we will run similar code to above but this time we will calculate the log fold change of TPM to see if there will be any differences with the data

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

# We see that there is no difference between the direction (TPM) and log fold change tables and that the trends seen in one, are seen in the other.

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
Viewable_trend_bar_plot <- ggplot(trend_counts, aes(x = Trend, y = n, fill = Species)) +
  geom_col(position = "dodge") +
  labs(
    title = "Gene Regulation Trends by Species (LFC 0h -> 96h)",
    x = "Regulation Trend",
    y = "Number of Genes"
  ) +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  theme(text = element_text(size = 12))

print(Viewable_trend_bar_plot)
```

![](Analysis_of_Human_and_Mouse_Gene_Expression_files/figure-gfm/Plotting%20Bar%20Plot-1.png)<!-- -->
\# We see from the bar plot that from the 0 time point to the 96 hour
time, there is more genes downregulated in Mouse and more genes
upregulated in Human after being treated with Doxycycline.
