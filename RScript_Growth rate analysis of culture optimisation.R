# ==================================================================================
# Growth rate analysis using growth curver
# ==================================================================================

# __________________________________________________________________________________
# Loading Packages 
# Data file: Dataset.RR.grwoth.curve.xlsx and named as growth_data
#___________________________________________________________________________________

library(readxl)
library(growthcurver)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(RColorBrewer)
library(broom)
library(multcompView)

growth_data <- growth_data %>%
  mutate( Condition = as.character(Condition), Factor = as.numeric(Factor), Time = as.numeric(Time),
        R1 = as.numeric(R1), R2 = as.numeric(R2), R3 = as.numeric(R3), R4 = as.numeric(R4) )

# Convert replicate columns to long format

data_long <- growth_data %>%
  pivot_longer( cols = c(R1, R2, R3, R4), names_to = "Replicate", values_to = "Value" ) %>%
  mutate( Replicate = factor(Replicate, levels = c("R1", "R2", "R3", "R4")) ) %>%
  filter(!is.na(Time),  !is.na(Value), !is.na(Condition), !is.na(Factor), is.finite(Time), is.finite(Value) ) %>%
  arrange( Condition, Factor, Replicate, Time )

# Summarised dataset of mean and SE across R1-R4

summary_data <- data_long %>%
  group_by(Condition, Factor,  Time) %>%
  summarise( n = sum(!is.na(Value)), Mean = mean(Value, na.rm = TRUE),
             SD = sd(Value, na.rm = TRUE), SE = SD / sqrt(n), .groups = "drop" )

# _________________________________________________________________________________
# Fit Growthcurver model to mean curve for every
# _________________________________________________________________________________

growth_fits <- data_long  %>%
  group_by(Condition, Factor, Replicate ) %>%
  nest() %>%
  mutate(fit = map(data,
      function(df) {
        
        df <- df %>%
          filter( !is.na(Time), !is.na(Value), is.finite(Time), is.finite(Value) ) %>%
          arrange(Time)
        
        # Need enough points for curve fitting
        if (nrow(df) < 5) {
          return(NULL)
        }
        
        tryCatch(
          SummarizeGrowth(  data_t = df$Time,  data_n = df$Value ),
          error = function(e) {
            message( "Growth model failed: ", conditionMessage(e)  )
   
            return(NULL)
          }
        )
      }
    )
  )

# Extract Growthcurver parameters

parameters_required <- c( "t_gen", "k", "k_se", "k_p", "n0", "n0_se", "n0_p",
  "r", "r_se", "r_p", "sigma", "df", "t_mid", "dt", "auc_l", "auc_e", "note")

growth_parameters <- growth_fits %>%
  mutate(parameters = map(  fit,
      function(gc) {
        if (is.null(gc)) {
          
          output <- as.list( rep( NA, length(parameters_required) ) )
          
          names(output) <- parameters_required
          output$note <- "Model not fitted"
          return(as_tibble(output))
        }
        
        vals <- as.list(gc$vals)
        missing_parameters <- setdiff( parameters_required, names(vals) )
        
        for (parameter in missing_parameters) {
          vals[[parameter]] <- NA
        }
        
        as_tibble(
          vals[parameters_required]
        )
      }
    )
  ) %>%
  select( Condition, Factor, Replicate, parameters ) %>%
  unnest(parameters) %>%
  mutate(  mu_max = log(2) / t_gen)

# Save parameters
write.csv(growth_parameters, "C:/Users/s11jb5/Desktop/Growthcurver_parameters.csv",
  row.names = FALSE)

# ______________________________________________________________________
# Plotting
# ______________________________________________________________________

#Palette

condition_palettes <- c(
  "Temperature"        = "Paired",
  "pH"                 = "Paired",
  "Ammonium_mM"        = "Dark2",
  "Urea_mM"            = "Set2",
  "Glucose_mM"         = "Set1",
  "Cellobiose_mM"      = "Set1",
  "Cellulose_percent"  = "Set1")

# Plotting function

plot_condition <- function(condition_name) {
  
  plot_data <- summary_data %>%
    filter(Condition == condition_name) %>%
    mutate(Factor = factor(Factor)) %>%
    arrange(Factor, Time)
  palette_name <- condition_palettes[[condition_name]]
  
  if (is.null(palette_name)) {
    stop(paste("No palette assigned for:", condition_name))
  }
  
  n_factors <- nlevels(plot_data$Factor)
  max_cols <- brewer.pal.info[palette_name, "maxcolors"]
  base_cols <- brewer.pal(min(max_cols, max(3, n_factors)),  palette_name)
  colours <- colorRampPalette(base_cols)(n_factors)
  
  ggplot( plot_data, aes( x = Time, y = Mean, colour = Factor, fill = Factor, group = Factor) ) +
    geom_ribbon( aes( ymin = Mean - SE, ymax = Mean + SE), alpha = 0.15,  colour = NA ) +
    geom_line( linewidth = 1 ) +
    geom_point( size = 2 ) +
    scale_colour_manual(values = colours) +
    scale_fill_manual(values = colours) +
    labs( title = condition_name, x = "Time", y = "Growth", colour = "Factor", fill = "Factor" ) +
    theme_classic(base_size = 14) +
    theme( plot.title = element_text( face = "bold", size = 16, hjust = 0.5),
      axis.title = element_text(face = "bold" ),
      axis.text = element_text( colour = "black" ),
      legend.title = element_text( face = "bold"),
      legend.position = "right" )
}


# Seperating plots

Temperature_plot <- plot_condition("Temperature")
pH_plot <- plot_condition("pH")
Ammonium_plot <- plot_condition("Ammonium_mM")
Urea_plot <- plot_condition("Urea_mM")
Glucose_plot <- plot_condition("Glucose_mM")
Cellobiose_plot <- plot_condition("Cellobiose_mM")
Cellulose_plot <- plot_condition("Cellulose_percent")

Temperature_plot
pH_plot
Ammonium_plot
Urea_plot
Glucose_plot
Cellobiose_plot
Cellulose_plot


# =======================================================================================
# Analyses of differences in predicted Growth rates from growthcurver and plotting
# =======================================================================================

# using growth_parameters as dataset

growth_parameters <- growth_parameters %>%
  mutate( Condition = as.character(Condition), Factor = as.numeric(Factor), r_day = as.numeric(r_day) ) %>%
  filter( !is.na(Condition), !is.na(Factor), !is.na(r_day), is.finite(Factor), is.finite(r_day))

# Summarised dataset of mean and SE across R1-R4

summary_rday <- growth_parameters %>%
  group_by(Condition, Factor) %>%
  summarise( n = n(), Mean = mean(r_day, na.rm = TRUE), SD = sd(r_day, na.rm = TRUE), SE = SD / sqrt(n), .groups = "drop" )

# ___________________________________________________________________________________
# Plotting
# ___________________________________________________________________________________

plot_theme <- theme_classic(base_size = 14) +
  theme(plot.title = element_text( face = "bold", size = 16, hjust = 0.5 ),
        axis.title = element_text( face = "bold"), axis.text = element_text( colour = "black") )

#Plot of bar plots

plot_bar <- function(condition_name) {
   plot_data <- summary_rday %>%
    filter(Condition == condition_name)
  
  ggplot(plot_data, aes(x = factor(Factor),y = Mean ) ) +
    geom_col( fill = "grey70", colour = "black",width = 0.7 ) +
    geom_errorbar( aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.15, linewidth = 0.7 ) +
    labs( title = condition_name, x = "Factor", y = expression(r~(day^{-1}))) +
     plot_theme
}

#Plot of barplots with smoothened line curve

plot_bar_smooth <- function(condition_name) {
  
  plot_data <- summary_rday %>%
    filter(Condition == condition_name) %>%
    arrange(Factor)
  
  factor_spacing <- diff(sort(unique(plot_data$Factor)))
  bar_width <- min(factor_spacing) * 0.25
  
  ggplot(plot_data, aes( x = Factor, y = Mean ) ) +
   # geom_line(colour = "grey30",linewidth = 0.7) +
    geom_col( fill = "grey90", colour = "black", width = bar_width, linewidth = 0.5) +
    geom_errorbar(aes( ymin = Mean - SE, ymax = Mean + SE), width = bar_width * 0.45, linewidth = 0.6 ) +
    geom_point( size = 2.5, shape = 21, fill = "grey50", colour = "black", stroke = 0.7) +
    labs( title = condition_name, x = "Factor", y = expression(r~(day^{-1})) ) +
    plot_theme
}

#Plot of dot plot with smoothened line curve

plot_dot_smooth <- function(condition_name) {
  
  plot_data <- summary_rday %>%
   filter(Condition == condition_name) %>%
    arrange(Factor)
  
  ggplot( plot_data, aes( x = Factor, y = Mean )) +
    geom_line(colour = "grey30",linewidth = 0.7) +
    geom_errorbar( aes( ymin = Mean - SE, ymax = Mean + SE ), width = 0.15, linewidth = 0.7, colour = "grey10" ) +
    geom_point( size = 3.5, shape = 21, fill = "grey70", colour = "black", stroke = 0.8 ) +
  # geom_smooth( method = "loess", formula = y ~ x, se = FALSE,  colour = "grey30", linewidth = 0.7, span = 0.8 ) +
    labs( title = condition_name, x = "Factor", y = expression(r~(day^{-1})) ) +
    
    plot_theme
}

# Seperating plots by each condition

conditions <- unique(summary_rday$Condition)
bar_plots <- setNames( lapply(conditions, plot_bar), conditions)
bar_smooth_plots <- setNames(lapply(conditions, plot_bar_smooth), conditions)
dot_smooth_plots <- setNames( lapply(conditions, plot_dot_smooth), conditions)

# Temperature
bar_plots$Temperature
bar_smooth_plots$Temperature
dot_smooth_plots$Temperature

# pH
bar_plots$pH
bar_smooth_plots$pH
dot_smooth_plots$pH

# Ammonium
bar_plots$Ammonium_mM
bar_smooth_plots$Ammonium_mM
dot_smooth_plots$Ammonium_mM

# Urea
bar_plots$Urea_mM
bar_smooth_plots$Urea_mM
dot_smooth_plots$Urea_mM

# Glucose
bar_plots$Glucose_mM
bar_smooth_plots$Glucose_mM
dot_smooth_plots$Glucose_mM

# Cellobiose
bar_plots$Cellobiose_mM
bar_smooth_plots$Cellobiose_mM
dot_smooth_plots$Cellobiose_mM

# Cellulose
bar_plots$Cellulose_percent
bar_smooth_plots$Cellulose_percent
dot_smooth_plots$Cellulose_percent


# ___________________________________________________________________________________
# Statistical analyses of differences in growth rates
# ___________________________________________________________________________________

# KRUSKAL-WALLIS test

stats_data <- growth_parameters %>%
  mutate(Condition = as.character(Condition), Factor_cat = factor(Factor) ) %>%
  filter( !is.na(Condition), !is.na(Factor_cat), !is.na(r_day) )

kruskal_results <- stats_data %>%
  group_by(Condition) %>%
  group_modify(~ { test <- kruskal.test( r_day ~ Factor_cat, data = .x )
    
 tibble( Chi_square = unname(test$statistic), df = unname(test$parameter), 
         p_value = test$p.value ) }) %>%
  ungroup() %>%
  mutate( significance = case_when( p_value < 0.001 ~ "***", p_value < 0.01  ~ "**", p_value < 0.05  ~ "*", 
      TRUE            ~ "ns" )  )

# DUNN POST-HOC 

dunn_results <- stats_data %>% group_by(Condition) %>%
  group_modify(~ {
    dunn <- FSA::dunnTest(  r_day ~ Factor_cat,  data = .x, method = "bh" )
    
    as_tibble(dunn$res) %>% transmute( Comparison = Comparison, Z = Z,
        P_unadjusted = P.unadj, P_adjusted = P.adj ) }) %>%
  ungroup() %>%
  mutate( significance = case_when( P_adjusted < 0.001 ~ "***", P_adjusted < 0.01  ~ "**",  P_adjusted < 0.05  ~ "*",
      TRUE               ~ "ns"  ) )

