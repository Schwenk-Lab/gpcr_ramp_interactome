#                     GPCR tissue expression analysis

#Load R packages:
library(tidyverse)
library(RColorBrewer)
library(scales)
library(patchwork)
library(rstatix)
library(pheatmap)
library(VennDiagram)

## Load expression data from https://www.proteinatlas.org/ (Version 23.0)
## Use the consensus data:
rna_tissue_data <- read.delim("rna_tissue_consensus.tsv", sep = "\t")
single_cell_data <- read.delim("rna_single_cell_type.tsv", sep = "\t")

# Read in table with information regarding the evidence levels for GPCR-RAMP interactions
file_name <- "2024-01-12_154035_cxdet_yes_no.csv"
seperator <- ","
protein_names_big <- read.csv(file_name, sep=seperator, dec=",", skip=0, stringsAsFactors=F, header=T)

# Extract the GPCR names
protein_names_big_names_only <- c(unique(protein_names_big$which_gpcr), "RAMP1", "RAMP2", "RAMP3")
## Change name of GPR1 to CMKLR2 as it is named in RNA_tissue_data
protein_names_big_names_only[grep(paste0("^GPR1$"), protein_names_big_names_only)] <- "CMKLR2"

# Extract the tissue expression data for the GPCRs in our data set:
rna_tissue_data_selected <- rna_tissue_data %>%
  filter(Gene.name %in% protein_names_big_names_only)

#------------------------------------- VENN DIAGRAMS  -----------------------------------------------------------------------------------
### Venn diagram for overlapping expected interactions for the different RAMPS:

ramp1 <- unname(c(protein_names_big %>%
                    filter(which_ramp == "RAMP1", class == "Yes") %>%
                    select(which_gpcr)))
ramp2 <- unname(c(protein_names_big %>%
                    filter(which_ramp == "RAMP2", class == "Yes") %>%
                    select(which_gpcr)))
ramp3 <- unname(c(protein_names_big %>%
                    filter(which_ramp == "RAMP3", class == "Yes") %>%
                    select(which_gpcr)))

venn_list <- list( "RAMP1" = unlist(ramp1),
                   "RAMP2" = unlist(ramp2),
                   "RAMP3" = unlist(ramp3))

# Create a Venn diagram from the list to pdf
myCol <- c("#BFEFFF", "#B1E063", "#EEAD0E")

# Create a Venn diagram from the list
venn_plot <- venn.diagram(
  x = venn_list,
  category.names = c("RAMP1", "RAMP2", "RAMP3"),
  filename = NULL,

  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

# Plot the Venn diagram
pdf("venn_plot_yes_interaction_test_test_test.pdf", width = 2, height = 2)
grid.draw(venn_plot)
dev.off()

#####                           Where we expect no interaction:

ramp1 <- unname(c(protein_names_big %>%
                    filter(which_ramp == "RAMP1", class == "No") %>%
                    select(which_gpcr)))
ramp2 <- unname(c(protein_names_big %>%
                    filter(which_ramp == "RAMP2", class == "No") %>%
                    select(which_gpcr)))
ramp3 <- unname(c(protein_names_big %>%
                    filter(which_ramp == "RAMP3", class == "No") %>%
                    select(which_gpcr)))

venn_list <- list( "RAMP1" = unlist(ramp1),
                   "RAMP2" = unlist(ramp2),
                   "RAMP3" = unlist(ramp3))


# Create a Venn diagram from the list
venn_plot <- venn.diagram(
  x = venn_list,
  category.names = c("RAMP1", "RAMP2", "RAMP3"),
  filename = NULL,
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

pdf("venn_plot_no_interaction.pdf", width = 2, height = 2)
grid.draw(venn_plot)
dev.off()


#------------------------------------- Scatter plots for number of GPCRs, and expected interaction percentage per RAMP3 ----------------------------------------------------

ramps <- c("RAMP1", "RAMP2", "RAMP3")

plot_list <- list()
result_summary_big <- NA
count_summary_big <- NA
for (r in ramps){
  
  protein_names_big_r <- protein_names_big %>%
    filter(which_ramp == r)
  
  ## Take only the cells where RAMP expression is above 1:
  cells_of_interest_r <- single_cell_data %>%
    filter(Gene.name %in% r) %>%
    filter(nTPM > 1) %>%
    select(Cell.type) %>%
    unlist()
  
  df_sc_r <- single_cell_data %>%
    filter(Gene.name %in% protein_names_big_names_only) %>% # Take out the GPCRs and the RAMPs
    filter(Cell.type %in% cells_of_interest_r) # Only for the cell types where we have the RAMP
  
  # Count the number of unique Gene.name for each Cell.type (# of GPCRs with expression above 1)
  gene_count <- df_sc_r %>%
    filter(nTPM > 1) %>% # GPCR must have expression > 1
    filter(Gene.name %in% protein_names_big_r$which_gpcr) %>% # Take out the GPCRs that interact or not interact with our specific RAMP (r)
    group_by(Cell.type) %>%
    summarise(total_genes = n_distinct(Gene.name))
  
  # Count the number of unique Gene.name for each Cell.type and class
  class_count <- df_sc_r %>%
    filter(nTPM > 1) %>%
    inner_join(protein_names_big_r, by = c("Gene.name" = "which_gpcr")) %>%
    group_by(Cell.type, class) %>%
    summarise(interaction_class = n_distinct(Gene.name))
  
  # Merge the counts
  df_r <- gene_count %>%
    left_join(class_count, by = "Cell.type") %>%
    mutate(
      percentage = interaction_class / total_genes * 100,
    ) 
  
  
  plot_2 <- ggplot(df_r, aes( x = reorder(Cell.type, -total_genes), y = total_genes)) +
    geom_point(aes(color = "# GPCRs"), position = position_dodge(width = 0.5), size = 3) +
    geom_point(aes(y = percentage, color = class), position = position_dodge(width = 0.5), size = 3) +
    labs(
      x = "Cell Type",
      y = "Total # GPCRs",
      title = paste0(r),
      color = ""
    ) +
    theme_classic() +
    scale_color_manual( values = c("#7570b3", "#d95f02", "#1b9e77")) +
    theme(axis.line = element_line(size = 1, colour = "black"), axis.text = element_text(size = 12), axis.title = element_text(size = 16), plot.title = element_text(size = 18), legend.title = element_text(size = 16), legend.text = element_text(size = 10), axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) +
    scale_y_continuous(sec.axis = sec_axis(~ ., name = "Percentage"))
  
  ### Calculate how much the percantage change per group:
  result_summary <- df_r %>%
    group_by(class) %>%
    summarise(
      max_percentage = max(percentage),
      min_percentage = min(percentage)
    )
  
  result_summary$ramp <- r
  result_summary_big <- rbind(result_summary_big, result_summary)
  
  
  result_summary <- df_r %>%
    summarise(
      total_genes_min = max(total_genes),
      total_genes_max = min(total_genes)
    )
  result_summary$ramp <- r
  count_summary_big <- rbind(count_summary_big, result_summary)
  
  plot_list[[r]] <- plot_2
  
}

plot_list[["RAMP1"]] / plot_list[["RAMP2"]] / plot_list[["RAMP3"]]


#------------------------------------- Violin plots for gpcr-ramp expression ratios in human tissues and cells ----------------------------------------------------

####                                     Prepare a table containing expression ratio for GPCR/RAMP in tissues:
ramps <- c("RAMP1", "RAMP2", "RAMP3")
tissues <- unique(rna_tissue_data_selected$Tissue)

ratio_result_table_tissue <- NA
for ( t in tissues){
  
  for ( ramp in ramps){
    proteins <- protein_names_big %>%
      filter(which_ramp == ramp)
    prot_names <- proteins$which_gpcr
    prot_names <- c(prot_names, ramp)
    
    ### Extract interaction_evidence:
    protein_names_ramp <- protein_names_big %>%
      filter(which_ramp == ramp )
    
    colnames(protein_names_ramp)[1] <- "Gene.name"
    
    # Extract rows where "Gene_name" is either the gpcr or ramp and calculate expression ratios
    filtered_data <- rna_tissue_data_selected %>%
      filter( Tissue == t  ) %>%
      filter(Gene.name %in% c( prot_names , ramp)) %>%
      filter( nTPM > 1) %>%
      mutate(nTPM_ratio = if (ramp %in% Gene.name) {
        nTPM / nTPM[Gene.name == ramp]
      } else {
        NA
      })
    
    filtered_data_evidence <- left_join(filtered_data, protein_names_ramp, by = "Gene.name")
    
    #Remove RAMP
    filtered_data_evidence_no_ramp <- filtered_data_evidence %>%
      filter( Gene.name != ramp) %>%
      mutate( interaction_evidence = case_when(class == "Yes" ~ "Yes",
                                               class == "No" ~ "No"))
    
    ratio_result_table_tissue <- rbind(ratio_result_table_tissue, filtered_data_evidence_no_ramp )
  } # end for ramp
} #end tissues

## Position for writing out the number of observations in the plot:
give.n <- function(x){
  return(c(y = -10, label = length(x))) 
}

boxplot_violin_t <- ggplot(ratio_result_table_tissue[2:nrow(ratio_result_table_tissue),], aes(x = which_ramp, y = log2(nTPM_ratio), fill = factor(class, levels = c("Yes", "No")), col = factor(class, levels = c("Yes", "No")) )) + # first row is NA in ratio_result_table
  geom_violin(alpha = 0.4) +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), col = "black", alpha = 0.4) +
  labs(
    title = "Expression ratio between GPCR and RAMP",
    subtitle = "Only keeping Tissues with RAMP and GPCR expression nTPM > 1",
    x = "",
    y = "log2(nTPM ratio (GPCR/RAMP))",
    fill = "Interaction Evidence",
    legend = "Interaction Evidence"
  ) +
  theme_classic() +
  scale_fill_brewer("Interaction Evidence" ,palette = "Dark2") +
  scale_color_brewer("Interaction Evidence",palette = "Dark2") + 
  stat_compare_means(label = "p.signif", method = "wilcox.test") +
  stat_summary(fun.data = give.n, geom = "text", fun = median, position = position_dodge(0.75)) +
  theme(panel.border = element_rect(color = "black", fill = NA),  text = element_text(family = "Helvetica"))


####                                           Prepare a table containing expression ratio for GPCR/RAMP in cell types:

ramps <- c("RAMP1", "RAMP2", "RAMP3")

## Only take out the cell types that have the ramps in them:
sc_data_selected <- single_cell_data %>%
  filter(Gene.name %in% ramps)
cells <- unique(sc_data_selected$Cell.type)

ratio_result_table <- NA
for (t in cells) {
  
  for (ramp in ramps) {
    
    proteins <- protein_names_big %>%
      filter(which_ramp == ramp)
    prot_names <- proteins$which_gpcr
    prot_names <- c(prot_names, ramp)
    
    ### Extract interaction_evidence:
    protein_names_ramp <- protein_names_big %>%
      filter(which_ramp == ramp )
    colnames(protein_names_ramp)[1] <- "Gene.name"
    
    # Extract rows where "Gene_name" is either the gpcrs or ramp and calculate expression ratios
    filtered_data <- single_cell_data %>%
      filter( Cell.type == t  ) %>%
      filter(Gene.name %in% c( prot_names , ramp)) %>%
      filter( nTPM > 1) %>%
      mutate(nTPM_ratio = if (ramp %in% Gene.name) {
        nTPM / nTPM[Gene.name == ramp]
      } else {
        NA
      })
    
    filtered_data_evidence <- left_join(filtered_data, protein_names_ramp, by = "Gene.name")
    
    #Remove RAMP
    filtered_data_evidence_no_ramp <- filtered_data_evidence %>%
      filter( Gene.name != ramp)
    
    ratio_result_table <- rbind(ratio_result_table, filtered_data_evidence_no_ramp )
  } # end for ramp
} #end cells


give.n <- function(x){
  return(c(y = -10, label = length(x))) 
}

boxplot_violin_sc <- ggplot(ratio_result_table[2:nrow(ratio_result_table),], aes(x =which_ramp, y = log2(nTPM_ratio), fill = factor(class, levels = c("Yes", "No")), col = factor(class, levels = c("Yes", "No")) )) + # first row is NA in ratio_result_table
  geom_violin(alpha = 0.4) +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), col = "black", alpha = 0.4) +
  labs(
    title = "Expression ratio between GPCR and RAMP",
    subtitle = "Cell types with RAMP and GPCR expression nTPM > 1",
    x = "",
    y = "log2(nTPM ratio (GPCR/RAMP))",
    fill = "Interaction Evidence",
    legend = "Interaction Evidence"
  ) +
  theme_classic() +
  scale_fill_brewer("Interaction Evidence" ,palette = "Dark2") +
  scale_color_brewer("Interaction Evidence",palette = "Dark2") + 
  stat_compare_means(label = "p.signif", method = "wilcox.test") +
  stat_summary(fun.data = give.n, geom = "text", fun.y = median, position = position_dodge(0.75)) +
  theme(panel.border = element_rect(color = "black", fill = NA),  text = element_text(family = "Helvetica"))


#------------------------------------- Heatmaps for ramp expression in different tissues, cell types  ----------------------------------------------------

# Binary heatmaps:

##  Tissue data, Use nTPM = 1 as cutoff for for expressed in tissue
ramp_t_d <- rna_tissue_data %>%
  filter(Gene.name %in% c("RAMP1", "RAMP2", "RAMP3")) %>%
  mutate(binary = ifelse(nTPM > 1, 1, 0))

heatmap_plot_tissue_t <- ramp_t_d %>%
  ggplot(aes(x = Gene.name, y = fct_rev(Tissue), fill = as.factor(binary))) +
  geom_tile(color = "white") +
  #scale_fill_manual(values = c("0" = "#EDD9A3", "1" = "#B1339E")) +  # Adjust colors as needed
  scale_fill_manual(values = c("0" = "#481567FF", "1" = "#3CBB75FF")) +  # Adju
  labs(
    title = "Tissue",
    x = "",
    y = "",
    fill = "Expression > 1"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  text = element_text(family = "Helvetica")) 


## Single cell, Use nTPM = 1 as cutoff for for expressed in cell type
ramp_sc_d <- single_cell_data %>%
  filter(Gene.name %in% c("RAMP1", "RAMP2", "RAMP3")) %>%
  mutate(binary = ifelse(nTPM > 1, 1, 0))

library(forcats)
heatmap_plot_sc_t <- ramp_sc_d %>%
  ggplot(aes(x = Gene.name, y = fct_rev(Cell.type), fill = as.factor(binary))) +
  geom_tile(color = "white") +
  #scale_fill_manual(values = c("0" = "#EDD9A3", "1" = "#B1339E")) +  # Adjust colors as needed
  scale_fill_manual(values = c("0" = "#481567FF", "1" = "#3CBB75FF")) +  # Adju
  labs(
    title = "Single cell",
    x = "",
    y = "",
    fill = "Expression > 1"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  text = element_text(family = "Helvetica")) 


### Heatmap for single-cell data log2 expression

## Reshape to wide format:
ramp_sc_d <- single_cell_data %>%
  filter(Gene.name %in% c("RAMP1", "RAMP2", "RAMP3"))

wide_data_ramp_sc <- reshape(ramp_sc_d[,c(2:4)], direction = "wide", 
                             idvar = "Gene.name", 
                             timevar = "Cell.type")

rownames(wide_data_ramp_sc) <- wide_data_ramp_sc$Gene.name
wide_data_ramp_sc[2:ncol(wide_data_ramp_sc)] <- apply(wide_data_ramp_sc[2:ncol(wide_data_ramp_sc)], 2, as.numeric)
rownames(wide_data_ramp_sc) <- wide_data_ramp_sc$Gene.name

# Log2-transform the expression data:
plot_df <- wide_data_ramp_sc %>%
  mutate_at(vars(-1), funs(ifelse(. > 1, log2(.), 0)))

## Remove nTPM from colname
colnames(plot_df) <- substr(colnames(plot_df), start = 6, stop = 100)

row_order <- c("RAMP1", "RAMP2", "RAMP3")

set.seed(1234)
d_ramp <- t(plot_df[, 2:ncol(plot_df)])
cc <- hclust(dist(t(d_ramp)))
ccdend <- as.dendrogram(cc) %>% reorder(c(2, 1, 2))
cc_new <- as.hclust(ccdend)
hm <- pheatmap(d_ramp, cluster_cols = cc_new, color = viridis::viridis(100)) 