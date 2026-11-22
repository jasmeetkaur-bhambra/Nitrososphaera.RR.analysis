#_______________________________________________________________________________
# PROTEOMICS BAR PLOTS 
# Plot of bars only when: BC vs C  > 0.5 & BC vs CC > 0.5
# B = biofilm 
# BC = B + C = biofilm with cellobiose
# C = FL = freeliving culture
# CC = FL + C = freeliving with cellobiose
# BC vs B
# CC vs C
# BC vs CC
#_______________________________________________________________________________

# Load dataset and packages
# Dataset = df = All significant proteins compiled

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(grid)

# FILTER PROTEINS: BC vs C  > 0.5 & BC vs CC > 0.5
# BC vs B does NOT determine selection.

df_selected <- df %>%
  filter( !is.na(logFC_BC_vs_C), !is.na(logFC_BC_vs_CC), logFC_BC_vs_C > 0.5,
    logFC_BC_vs_CC > 0.5, !is.na(Final_annotation), Final_annotation != "",
   !is.na(Category), Category != "", !is.na(Major_Category), Major_Category != "" ) %>%
   mutate( Protein_ID = row_number() )

# Convert data to long format for plotting removing NAs

df_long <- df_selected %>% 
  select(Accession, Protein_ID, Major_Category, Category, Final_annotation, logFC_BC_vs_C,
          logFC_BC_vs_CC, logFC_BC_vs_B) %>%
  pivot_longer(cols = c( logFC_BC_vs_C, logFC_BC_vs_CC,logFC_BC_vs_B),
    names_to = "Comparison", values_to = "logFC") %>%
  mutate(Comparison = recode( Comparison, "logFC_BC_vs_C"  = "BC vs C",
      "logFC_BC_vs_CC" = "BC vs CC", "logFC_BC_vs_B"  = "BC vs B") ) %>%
  filter( !is.na(logFC))

# Set order & colours

df_long <- df_long %>%
  mutate( Comparison = factor( Comparison, levels = c("BC vs B", "BC vs CC",
        "BC vs C" ) ) )

comparison_colours <- c("BC vs C"  = "#a2dff7", "BC vs CC" = "#5a71af",
  "BC vs B"  = "#021f26")


# DEFINE MAJOR CATEGORY GROUPS & ORDER

plot1_categories <- c( "Protein Quality Control", "Defense & Stress Response",
  "Stress Protection & Homeostasis", "Energy & Respiratory Metabolism",
  "Cell Division & Growth")

plot2a_categories <- c( "Biofilm Formation")

plot2b_categories <- c("Cell Envelope & Membrane Adaptation", "Carbohydrate metabolism")

plot3_categories <- c( "Coenzyme biosynthesis & metabolism", "Amino Acid Metabolism")

plot4_categories <- c("Information Processing & Gene Regulation")

plot5_categories <- c( "Transport", "Environmental Sensing & Regulation")

plot6_categories <- c( "Uncharacterised")

major_category_order <- c( plot1_categories, plot2a_categories, plot2b_categories,
  plot3_categories, plot4_categories, plot5_categories, plot6_categories)

df_long <- df_long %>%
    mutate(Major_Category = factor( Major_Category, levels = major_category_order))


# PREPARE DATA FOR PLOTTING

prepare_plot_data <- function(data) {
  
  protein_order <- data %>%
     filter(Comparison == "BC vs C" ) %>%
    arrange( Major_Category, Category, logFC ) %>%
    pull(Protein_ID) %>%
    unique()
  
  data <- data %>%
     mutate( Gene_plot = paste0(Final_annotation,  "___", Protein_ID ))
  
  gene_levels <- data %>% 
    distinct( Protein_ID, Gene_plot ) %>%
    mutate( Protein_order = match( Protein_ID, protein_order )) %>%
    arrange( Protein_order ) %>%
    pull( Gene_plot)

  data <- data %>%
     mutate( Gene_plot = factor( Gene_plot, levels = gene_levels ) )
  
  return(data)
}


# PLOTTING FUNCTION

make_protein_plot <- function( data, title_text, facet_ncol = 1
) {
  
  plot_data <- prepare_plot_data(data)
  
  p <- ggplot( plot_data,  aes( x = logFC, y = Gene_plot, fill = Comparison,
      group = Comparison ) ) +
  geom_col( position = position_dodge( width = 0.85 ), width = 0.78, ,
            colour = "grey10", linewidth = 0.25 ) +
  geom_vline( xintercept = 0.5, linetype = "dashed", colour = "grey70",
    linewidth = 0.4 ) +
    facet_grid( Category ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = comparison_colours, breaks = c("BC vs B",
      "BC vs CC", "BC vs C" ), drop = FALSE  ) +
  scale_y_discrete( labels = function(y) {
      sub( "___.*$", "", y )
       },
    drop = TRUE ) +
  labs( title = title_text,  x = "logFC", y = NULL, fill = "Comparison" ) +
   theme_classic( base_size = 11 ) +
  theme(  axis.text.y = element_text( size = 6,  hjust = 1 ),
        axis.text.x = element_text( size = 9 ),
      axis.title.x = element_text( size = 10),
      strip.text.y.right = element_text( size = 9, face = "bold", angle = 0),
     strip.background = element_rect( fill = "grey95", colour = "grey5", linewidth = 0.3 ),
      panel.spacing.y = unit( 1, "lines"), panel.spacing.x = unit( 1, "lines"),
      legend.position = "bottom", legend.title = element_text( size = 9 ),
     legend.text = element_text( size = 8 ),  plot.margin = margin( 10, 10,  10, 10)
     
    )
  
  return(p)
  
}


# PLOTTING

data_plot1 <- df_long %>% filter(Major_Category %in% plot1_categories ) %>%
  mutate( Major_Category = factor( Major_Category,   levels = plot1_categories )  )
make_protein_plot(data_plot1, "Protein quality, stress, energy and growth")

data_plot2a <- df_long %>% filter( Major_Category == "Biofilm Formation")
make_protein_plot( data_plot2a,  "Biofilm Formation")

data_plot2b <- df_long %>% filter(  Major_Category %in% plot2b_categories) %>%
    mutate(  Major_Category = factor(Major_Category,  levels = plot2b_categories  ) )
make_protein_plot( data_plot2b,  "Cell envelope and carbohydrate metabolism")

data_plot3 <- df_long %>% filter( Major_Category %in% plot3_categories ) %>%
   mutate( Major_Category = factor( Major_Category, levels = plot3_categories ))
make_protein_plot(data_plot3, "Coenzyme and amino acid metabolism")

data_plot4 <- df_long %>% filter( Major_Category %in% plot4_categories ) %>%
  mutate( Major_Category = factor( Major_Category, levels = plot4_categories ))
make_protein_plot(data_plot4, "Information Processing & Gene Regulation")

data_plot5 <- df_long %>% filter( Major_Category %in% plot5_categories ) %>%
  mutate( Major_Category = factor( Major_Category, levels = plot5_categories ))
make_protein_plot(data_plot5, "Transport and environmental responses")

# TABLE OF SELECTED PROTEINS & BARS

selected_proteins <- df_selected %>%
  select(Major_Category, Category, Final_annotation, logFC_BC_vs_C,
    logFC_BC_vs_CC, logFC_BC_vs_B ) %>%
  arrange(Major_Category, Category, Final_annotation)

bar_check <- df_long %>%
  count(Accession,Protein_ID, Final_annotation, name = "Number_of_bars")
table(bar_check$Number_of_bars)
