#_______________________________________________________________________________
# PATHWAY-LEVEL PIE BUBBLE PLOTS
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
library(scatterpie)

df_long <- df %>% select(Major_Category, logFC_BC_vs_B, logFC_CC_vs_C, 
          logFC_BC_vs_CC, logFC_pH5_vs_pH6, logFC_pH7_vs_pH6, logFC_pH5_vs_pH7) %>%
     pivot_longer(cols = c(logFC_BC_vs_B, logFC_CC_vs_C, logFC_BC_vs_CC,
               logFC_pH5_vs_pH6, logFC_pH7_vs_pH6, logFC_pH5_vs_pH7),
    names_to = "Comparison", values_to = "logFC") %>%
  mutate(Comparison = recode( Comparison, "logFC_BC_vs_B"  = "BC vs B",
      "logFC_CC_vs_C"  = "CC vs C", "logFC_BC_vs_CC" = "BC vs CC",  
      "logFC_pH5_vs_pH6" = "pH5 vs pH6", "logFC_pH7_vs_pH6" = "pH7 vs pH6",
      "logFC_pH5_vs_pH7" = "pH5 vs pH7"),
    FC_group = case_when(logFC > 0.5  ~ "Upregulated", logFC < -0.5 ~ "Downregulated")) %>%
  filter(!is.na(Major_Category), !is.na(logFC),
    !grepl("uncharacter", Major_Category, ignore.case = TRUE))

df_long$Comparison <- factor(df_long$Comparison, levels = c("BC vs B", "CC vs C",
                   "BC vs CC", "pH5 vs pH6", "pH7 vs pH6", "pH5 vs pH7"  ))

# Pie data

prepare_pie_data <- function(data) {
  
  summary_data <- data %>% count(Major_Category, FC_group) %>%
  complete(Major_Category, FC_group = c("Upregulated", "Downregulated"),
   fill = list(n = 0)) %>%
    pivot_wider(names_from = FC_group, values_from = n, values_fill = 0 ) %>%
   mutate( Total = Upregulated + Downregulated )
  
category_order <- summary_data %>% arrange(Total) %>% pull(Major_Category)
 summary_data <- summary_data %>%
    mutate( Major_Category = factor( Major_Category, levels = category_order),
        y = as.numeric(Major_Category), x = 1)
  
  return(summary_data)
}

pie_BC_vs_B <- df_long %>% filter(Comparison == "BC vs B") %>% prepare_pie_data()
pie_CC_vs_C <- df_long %>% filter(Comparison == "CC vs C") %>% prepare_pie_data()
pie_BC_vs_CC <- df_long %>% filter(Comparison == "BC vs CC") %>% prepare_pie_data()
pie_pH5_vs_pH6 <- df_long %>% filter(Comparison == "pH5 vs pH6") %>% prepare_pie_data()
pie_pH7_vs_pH6 <- df_long %>% filter(Comparison == "pH7 vs pH6") %>% prepare_pie_data()
pie_pH5_vs_pH7 <- df_long %>% filter(Comparison == "pH5 vs pH7") %>% prepare_pie_data()

plot_pie_bubble <- function(plot_data, title_text) {
 
 plot_data <- plot_data %>% mutate(radius = sqrt(Total) * 0.03)

 ggplot() +
  geom_scatterpie(data = plot_data, aes(x = x, y = y, r = radius),
    cols = c( "Upregulated", "Downregulated" ), colour = "black", linewidth = 0.25) +
  scale_fill_manual(values = c("Upregulated" = "seagreen", "Downregulated" = "pink" ),
   breaks = c( "Upregulated", "Downregulated"), labels = c( "logFC > 0.5", "logFC < -0.5"),
   name = "Protein response") +
   scale_x_continuous( breaks = 1, labels = title_text, limits = c(0.5, 1.5)) +
    scale_y_continuous(breaks = seq_along( levels(plot_data$Major_Category)),
    labels = levels( plot_data$Major_Category ) ) +
    coord_fixed() +
    labs(title = title_text, x = NULL, y = NULL, fill = "logFC") +
     theme_classic(base_size = 12) +
    theme( axis.text.x = element_text( size = 11, face = "bold"),
        axis.text.y = element_text(size = 9), axis.ticks = element_blank(),
       legend.position = "right", plot.title = element_text( face = "bold",
        size = 13) )
}

p_BC_vs_B <- plot_pie_bubble( pie_BC_vs_B, "BC vs B")
p_BC_vs_B
p_CC_vs_C <- plot_pie_bubble( pie_CC_vs_C, "CC vs C")
p_CC_vs_C
p_BC_vs_CC <- plot_pie_bubble(pie_BC_vs_CC, "BC vs CC")
p_BC_vs_CC
p_pH5_vs_pH6 <- plot_pie_bubble( pH5_vs_pH6, "pH5 vs pH6")
p_pH5_vs_pH6
p_pH7_vs_pH6 <- plot_pie_bubble( pie_pH7_vs_pH6, "pH7 vs pH6")
p_pH7_vs_pH6
p_pH5_vs_pH7 <- plot_pie_bubble(pie_pH5_vs_pH7, "pH5 vs pH7")
p_pH5_vs_pH7
