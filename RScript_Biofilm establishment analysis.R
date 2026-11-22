# ______________________________________________________________________________
# BIOFILM CULTURE ANALYSIS
# ______________________________________________________________________________

# PACKAGES

library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(growthcurver)
library(openxlsx)
library(car)
library(emmeans)


# Data: Dataset.biofilm.cultivation and named db

#Ordering 

combo_order <- c("Transfer 0 Unamended", "Transfer 0 Cellobiose",
  "Transfer 1 Unamended", "Transfer 1 Cellobiose", "Transfer 2 Unamended",
  "Transfer 2 Cellobiose", "Transfer 3 Unamended", "Transfer 3 Cellobiose",
  "Transfer 4 Unamended", "Transfer 4 Cellobiose", "Transfer 5 Unamended",
  "Transfer 5 Cellobiose")

db <- db %>%
  mutate(Combo.condition = factor( Combo.condition, levels = combo_order),
         Transfer = factor(Transfer, levels = paste0("Transfer ", 0:5)),
         Condition = factor(Condition, levels = c("Unamended", "Cellobiose")))

# assigning colours

combo_colours <- c("Transfer 0 Unamended" = "#87b3d9",
  "Transfer 0 Cellobiose" = "gold4", "Transfer 1 Unamended" = "#3e8b96",
  "Transfer 1 Cellobiose" = "#feba4d", "Transfer 2 Unamended" = "#b3d7a5",
  "Transfer 2 Cellobiose" = "#f69bba", "Transfer 3 Unamended" = "#7ab58a",
  "Transfer 3 Cellobiose" = "#f13e8e", "Transfer 4 Unamended" = "#69ab56",
  "Transfer 4 Cellobiose" = "#a32a6b", "Transfer 5 Unamended" = "#0a4c25",
  "Transfer 5 Cellobiose" = "#7f1c26")

# Summarise dataset for mean and SE

summary_se <- function(data, time_variable) {
   data %>%
    group_by(Combo.condition, Transfer, Condition, {{ time_variable }}) %>%
    summarise(n = sum(!is.na(Nitrite)), mean_Nitrite = mean(Nitrite, na.rm = TRUE),
      sd_Nitrite = sd(Nitrite, na.rm = TRUE), se_Nitrite = sd_Nitrite / sqrt(n),
      .groups = "drop" )
}

summary_cont <- summary_se(db, Time.cont)
summary_time <- summary_se(db, Time)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Plotting the dataset
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Basic theme for plotting 

plot_theme <- theme_classic(base_size = 13) +
theme(axis.title = element_text(size = 14), legend.title = element_blank(),
   axis.text = element_text(size = 12), legend.text = element_text(size = 10),
   strip.text = element_text( size = 13, face = "bold"),
   panel.spacing = unit(1, "lines"))


# PLOT of Time.continuous vs Nitrite

ggplot(summary_cont, aes(x = Time.cont, y = mean_Nitrite,
    colour = Condition, group = Combo.condition)) +
  geom_errorbar(aes( ymin = mean_Nitrite - se_Nitrite,
      ymax = mean_Nitrite + se_Nitrite), width = 0.25, linewidth = 0.6) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  scale_colour_manual(values = fill_colours, drop = FALSE) +
  labs( x = "Continuous time", y = "Nitrite") + plot_theme


# PLOT of Time vs Nitrite

ggplot(summary_time, aes(x = Time, y = mean_Nitrite,
    colour = Combo.condition, group = Combo.condition)) +
  geom_errorbar(aes( ymin = mean_Nitrite - se_Nitrite,
      ymax = mean_Nitrite + se_Nitrite), width = 0.20, linewidth = 0.6) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  scale_colour_manual(values = combo_colours, drop = FALSE) +
  labs(x = "Time", y = "Nitrite") + plot_theme


# PLOT of Time vs Nitrite subset by condition

ggplot(summary_time, aes(x = Time, y = mean_Nitrite,
    colour = Combo.condition, group = Combo.condition)) +
  geom_errorbar(aes(ymin = mean_Nitrite - se_Nitrite,
      ymax = mean_Nitrite + se_Nitrite), width = 0.20, linewidth = 0.6) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  facet_wrap( ~ Condition, ncol = 2) +
  scale_colour_manual(values = combo_colours, drop = FALSE) +
  labs(x = "Time", y = "Nitrite") + plot_theme

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GROWTHCURVER Modelling 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Growth curves are fitted for each condition, transfer & replicate for Time

fit_growthcurve <- function(data) {
  d <- data %>% filter( !is.na(Time), !is.na(Nitrite)) %>% arrange(Time)
  # To remove curves with values less than 4
  if (nrow(d) < 4) {
    return(
      tibble(k = NA_real_, n0 = NA_real_, r = NA_real_, t_mid = NA_real_,
        t_gen = NA_real_, auc_l = NA_real_, auc_e = NA_real_, sigma = NA_real_,
        note = "Too few observations")
    )
  }
  
  fit <- tryCatch(
      SummarizeGrowth( data_t = d$Time, data_n = d$Nitrite),
    error = function(e) NULL)

  if (is.null(fit)) {
    return( tibble( k = NA_real_, n0 = NA_real_, r = NA_real_, t_mid = NA_real_,
        t_gen = NA_real_, auc_l = NA_real_, auc_e = NA_real_, sigma = NA_real_,
        note = "Fit failed"))
  }
  # Extracting metrics
  metrics <- fit$vals
  tibble(k = metrics$k, n0 = metrics$n0, r = metrics$r, t_mid = metrics$t_mid,
         t_gen = metrics$t_gen, auc_l = metrics$auc_l, auc_e = metrics$auc_e,
         sigma = metrics$sigma, note = "Fit successful" )
}

growth_parameters <- db %>%
  group_by(Transfer, Condition, Combo.condition, Replicate) %>%
  nest() %>% mutate( growth_fit = map( data, fit_growthcurve)) %>%
  select( -data) %>% unnest( growth_fit ) %>%  ungroup()


# Exporting to excel

wb <- createWorkbook()
addWorksheet( wb, "Growth parameters")
writeData( wb, "Growth parameters", growth_parameters)
saveWorkbook( wb, "Biofilm_GrowthCurver_Parameters.xlsx", overwrite = TRUE)


# PLOT ALL GROWTHCURVER FITS
# Create folder for individual plots
plot_folder <- "GrowthCurver_Fitted_Curves"
if (!dir.exists(plot_folder)) {
  dir.create(plot_folder)
}

#Function to fit growthcurver data

make_growth_fit <- function(data)   {
  
  d <- data %>% filter( !is.na(Time), !is.na(Nitrite) ) %>% arrange(Time)
  # Need enough observations
  if (nrow(d) < 4) {
    return(NULL)
  }
  # Fit Growthcurver
  fit <- tryCatch(
    SummarizeGrowth( data_t = d$Time, data_n = d$Nitrite),
    error = function(e) NULL)
  if (is.null(fit)) {
    return(NULL)
  }
  # Extract fitted parameters
  k  <- fit$vals$k
  n0 <- fit$vals$n0
  r  <- fit$vals$r
  # Fine time sequence for smooth fitted curve
  new_time <- seq(min(d$Time), max(d$Time), length.out = 300)
  # Logistic equation used by Growthcurver
  fitted_nitrite <- k / ( 1 + ((k - n0) / n0) * exp(-r * new_time)
  tibble( Time = new_time, Fitted_Nitrite = fitted_nitrite )
}

curve_data <- db %>%
  group_by( Transfer, Condition, Combo.condition, Replicate) %>%
  nest() %>%
  mutate(fit_data = map( data, make_growth_fit)) %>%
  ungroup()

all_growth_plots <- pmap(
  curve_data,
  function( Transfer, Condition, Combo.condition, Replicate, data, fit_data
  ) {
    make_growth_plot( raw_data = data, fitted_data = fit_data, transfer = Transfer,
      condition = Condition, combo = Combo.condition, replicate = Replicate )
  }
)

# 1 plot per PNG

for (i in seq_along(all_growth_plots)) {
  transfer_i <- as.character( curve_data$Transfer[i])
  condition_i <- as.character( curve_data$Condition[i])
  replicate_i <- as.character( curve_data$Replicate[i] )
  filename <- paste0( plot_folder, "/", gsub(" ", "_", transfer_i), "_",
    gsub(" ", "_", condition_i), "_Replicate_", replicate_i, ".png")
  
  ggsave( filename = filename, plot = all_growth_plots[[i]], width = 6,
    height = 5, dpi = 600 )
}


if (!dir.exists(plot_folder)) { dir.create(plot_folder, recursive = TRUE)}

for (i in seq_len(nrow(curve_data))) {
  transfer_i <- as.character(curve_data$Transfer[i])
  condition_i <- as.character(curve_data$Condition[i])
  combo_i <- as.character(curve_data$Combo.condition[i])
  replicate_i <- as.character(curve_data$Replicate[i])
  raw_i <- curve_data$data[[i]]
  fit_i <- curve_data$fit_data[[i]]
  p <- make_growth_plot(
    raw_data = raw_i, fitted_data = fit_i, transfer = transfer_i,
    condition = condition_i, combo = combo_i, replicate = replicate_i )

  filename <- file.path(plot_folder, paste0( gsub(" ", "_", transfer_i),
      "_", gsub(" ", "_", condition_i), "_Replicate_", replicate_i, ".png" ) )

  png(filename = filename, width = 1800, height = 1500, res = 300)
  print(p)
  dev.off()
  
  cat( "Saved:", transfer_i, "|", condition_i, "| Replicate", replicate_i, "\n")
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot and statistics of all extracted parameters
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# rewording dataset growth_parameters as df

df <- df %>%
  mutate( Transfer = factor( Transfer, levels = paste0("Transfer ", 0:5)),
    Condition = factor( Condition, levels = c( "Unamended", "Cellobiose" ) ) )

# Setting colours for plotting

fill_colours <- c( "Cellobiose" = "#b47cab", "Unamended"  = "#7cb486")
outline_colours <- c( "Cellobiose" = "#72426A", "Unamended"  = "#2d6c25")

# Calculating mean and SE per replicate per condition and transfer

summary_df <- df %>% filter(note == "Fit successful") %>%
  group_by( Transfer, Condition) %>%
  summarise( n_k = sum(!is.na(k)), mean_k = mean(k, na.rm = TRUE),
    se_k = sd(k, na.rm = TRUE) / sqrt(n_k),
    n_r = sum(!is.na(r)), mean_r = mean(r, na.rm = TRUE),
    se_r = sd(r, na.rm = TRUE) / sqrt(n_r),
    n_tmid = sum(!is.na(t_mid)), mean_tmid = mean(t_mid, na.rm = TRUE),
    se_tmid = sd(t_mid, na.rm = TRUE) / sqrt(n_tmid),
    n_tgen = sum(!is.na(t_gen)), mean_tgen = mean(t_gen, na.rm = TRUE),
    se_tgen = sd(t_gen, na.rm = TRUE) / sqrt(n_tgen),
    .groups = "drop" )

# Plot common theme

plot_theme <- theme_classic(base_size = 14) +
   theme(axis.title = element_text( size = 14),axis.text = element_text(size = 12),
   legend.title = element_blank(), legend.text = element_text( size = 12),
    legend.position = "top")


# PLOT of r / Growth rate

ggplot(summary_df,
  aes(y = mean_r,x = Transfer, fill = Condition, colour = Condition )) +
  geom_col( position = position_dodge(width = 0.7), width = 0.6, linewidth = 0.8) +
  geom_errorbar(aes( ymin = mean_r - se_r, ymax = mean_r + se_r),
    position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.8) +
    scale_fill_manual(values = fill_colours) + scale_colour_manual(values = outline_colours) +
  geom_point(data = df, aes(y = r,x = Transfer, fill = Condition, colour = Condition),
    inherit.aes = FALSE, size = 2, position = position_dodge(width = 0.7) ) +
  labs( y = "r", x = "" ) + plot_theme


# PLOT of t-mid / Time at midpoint of logistic growth

ggplot(summary_df,
       aes(x = mean_tmid,y = Transfer, fill = Condition, colour = Condition )) +
  geom_col( position = position_dodge(width = 0.7), width = 0.6, linewidth = 0.8) +
  geom_errorbar(aes( xmin = mean_tmid - se_tmid, xmax = mean_tmid + se_tmid),
                position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.8) +
  geom_point(data = df, aes(x = t_mid,y = Transfer, fill = Condition, colour = Condition),
             inherit.aes = FALSE, size = 2, position = position_dodge(width = 0.7) ) +
  scale_fill_manual(values = fill_colours) + scale_colour_manual(values = outline_colours) +
  scale_y_discrete(  limits = rev(levels(summary_df$Transfer)) ) +
  labs( x = "t-mid", y = "" ) + plot_theme


# PLOT of t-gen / Generation time

ggplot(summary_df,
       aes(x = mean_tgen,y = Transfer, fill = Condition, colour = Condition )) +
  geom_col( position = position_dodge(width = 0.7), width = 0.6, linewidth = 0.8) +
  geom_errorbar(aes( xmin = mean_tgen - se_tgen, xmax = mean_tgen + se_tgen),
                position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.8) +
  geom_point(data = df, aes(x = t_gen,y = Transfer, fill = Condition, colour = Condition),
             inherit.aes = FALSE, size = 2, position = position_dodge(width = 0.7) ) +
  scale_fill_manual(values = fill_colours) + scale_colour_manual(values = outline_colours) +
  scale_y_discrete(  limits = rev(levels(summary_df$Transfer)) ) +
  labs( x = "Generation time", y = "" ) + plot_theme


# STATISTICAL ANALYSIS OF EXTRACTED GROWTHCURVER PARAMETERS

options(contrasts = c("contr.sum", "contr.poly" ))

# Models

model_r <- lm( r ~ Transfer * Condition, data = df)
model_tmid <- lm( t_mid ~ Transfer * Condition, data = df)
model_tgen <- lm(t_gen ~ Transfer * Condition, data = df)

# Model fit

par(mfrow = c(2, 2))

plot(model_r)
plot(model_tmid)
plot(model_tgen)

# SHAPIRO-WILK

shapiro_r <- shapiro.test( residuals(model_r))
shapiro_tmid <- shapiro.test( residuals(model_tmid))
shapiro_tgen <- shapiro.test( residuals(model_tgen))

shapiro_r
shapiro_tmid
shapiro_tgen

# HOMOGENEITY OF VARIANCE

levene_r <- leveneTest(r ~ interaction(Transfer, Condition), data = df)
levene_tmid <- leveneTest(t_mid ~ interaction(Transfer, Condition), data = df)
levene_tgen <- leveneTest(t_gen ~ interaction(Transfer, Condition), data = df)

levene_r
levene_tmid
levene_tgen

# ANOVA & TUKEY test

anova_r <- aov(model_r)
summary(anova_r)
TukeyHSD(anova_r)

anova_tmid <- aov(model_tmid)
summary(anova_tmid)
TukeyHSD(anova_tmid)

anova_tgen <- aov(model_tgen)
summary(anova_tgen)
TukeyHSD (anova_tgen)

#
#
#

#_______________________________________________________________________________
# Plot and statistics of cell density 
#_______________________________________________________________________________

# Embedding dataset within script

dc <- data.frame( Condition = c( "Unamended", "Unamended", "Unamended",
    "Cellobiose", "Cellobiose", "Cellobiose"),
  Cells = c( 2.31E6, 7.60E5, 2.19E6, 4.27E6, 4.01E6, 5.37E6 ))

dc$Condition <- factor( dc$Condition, levels = c("Unamended", "Cellobiose"))

fill_colours <- c("Cellobiose" = "#b47cab", "Unamended"  = "#7cb486")
outline_colours <- c("Cellobiose" = "#72426A", "Unamended"  = "#2d6c25")


# MEan and SE summary

summary_dc <- dc %>% group_by(Condition) %>%
  summarise( mean_cells = mean(Cells), sd_cells   = sd(Cells),
    n = n(), se_cells   = sd_cells / sqrt(n),
    .groups = "drop" )

# T-test

t.test(Cells ~ Condition, data = dc)

# Plot

ggplot(summary_dc, aes( x = mean_cells, y = Condition, fill = Condition,
    colour = Condition)) + geom_col(width = 0.6, linewidth = 0.8 ) +
  geom_errorbar(aes(  xmin = mean_cells - se_cells, xmax = mean_cells + se_cells),
    width = 0.18, linewidth = 0.8 ) +
  geom_point( data = dc, aes( x = Cells, y = Condition, colour = Condition ),
    inherit.aes = FALSE, size = 2, position = position_jitter(height = 0.07)) +
  scale_fill_manual(values = fill_colours) + scale_colour_manual(values = outline_colours) +
  scale_x_continuous(labels = scales::label_scientific() ) +
  labs( x = expression("Cell density (cells " * mL^{-1} * ")"), y = NULL) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none", axis.text = element_text(colour = "black"), 
        axis.title = element_text(colour = "black"))

