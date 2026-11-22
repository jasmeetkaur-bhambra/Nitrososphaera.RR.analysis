#_______________________________________________________________________________
# BIOFILM pH GROWTH CURVES
#_______________________________________________________________________________

# Input: Dataset_Biofilm.pH.range.test & dataset named as df

#Loading packages 

library(tidyverse)
library(growthcurver)
library(tidyverse)
library(car)
library(emmeans)
library(dplyr)
library(FSA)

# Loading and correcting dataset

df <- df %>%
  mutate(Condition = as.factor(Condition), pH = as.factor(pH),
    Replicate = as.factor(Replicate), Condition.2 = as.factor(Condition.2),
    Time.d = as.numeric(Time.d), Nitrite.uM = as.numeric(Nitrite.uM) )


# Summarised dataset of mean and SE

summary_df <- df %>% group_by(Condition, Condition.2, Time.d) %>%
  summarise( n = sum(!is.na(Nitrite.uM)), mean_nitrite = mean(
      Nitrite.uM, na.rm = TRUE ),
     sd_nitrite = sd( Nitrite.uM, na.rm = TRUE ),
    se_nitrite = sd_nitrite / sqrt(n),
    .groups = "drop")

# Subsetting colours for plots

condition_colours <- c("cellobiose5"   = "#8f2930", "cellobiose5.5" = "#eb6669",
  "cellobiose6"   = "#f99bba","unamended6"    = "#95e0e7",
  "unamended5.5"  = "#51b8d5",  "unamended5"    = "#2a4d66")

df$Condition.2 <- factor(df$Condition.2,  levels = condition_order)


# Basic plot theme

plot_theme <- theme_classic(base_size = 14) +
  theme(axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text( size = 12, colour = "black" ),
    axis.line = element_line( colour = "black", linewidth = 0.7 ),
    axis.ticks = element_line( colour = "black" ),
    legend.title = element_blank(), legend.text = element_text(size = 11 ),
    plot.title = element_text( size = 14, face = "bold", hjust = 0.5 ))

# PLOT: mean nitrite over time in 1 single plot

ggplot( summary_df, aes( x = Time.d, y = mean_nitrite, colour = Condition.2,
    group = Condition.2 )) + geom_line( linewidth = 1 ) + geom_point( size = 3) +
  geom_errorbar( aes( ymin = mean_nitrite - se_nitrite,
         ymax = mean_nitrite + se_nitrite), width = 0.20, linewidth = 0.6 ) +
  scale_colour_manual( values = condition_colours, drop = FALSE ) +
  labs(x = "Time (days)", y = expression(Nitrite~(mu*M))) +  plot_theme

# PLOT: mean nitrite over time subset by condition

ggplot( summary_df, aes( x = Time.d, y = mean_nitrite, colour = Condition.2,
                         group = Condition.2 )) + geom_line( linewidth = 1 ) + geom_point( size = 3) +
  geom_errorbar( aes( ymin = mean_nitrite - se_nitrite,
                      ymax = mean_nitrite + se_nitrite), width = 0.20, linewidth = 0.6 ) +
  facet_wrap( ~ Condition, nrow = 1) +
  scale_colour_manual( values = condition_colours, drop = FALSE ) +
  labs(x = "Time (days)", y = expression(Nitrite~(mu*M))) +  plot_theme


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GROWTHCURVER analysis of growth
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fit_growthcurve <- function(data) {
  
  data <- data %>% filter( !is.na(Time.d), !is.na(Nitrite.uM) ) %>% arrange(Time.d)

  if(nrow(data) < 5) {
    return(NULL)
  }

  fit <- tryCatch( SummarizeGrowth( data$Time.d, data$Nitrite.uM ),
    error = function(e) NULL)
  
  return(fit)
}


# Splitting data by each replicate

split_df <- df %>% group_by(  Condition,  pH,  Replicate, Condition.2) %>%
  group_split()


# Fitting the curves

all_parameters <- list()

for(i in seq_along(split_df)) {
  
  dat <- split_df[[i]]
  
  condition <- as.character(unique(dat$Condition))
  ph <- as.character( unique(dat$pH))
  replicate <- as.character( unique(dat$Replicate))
  condition2 <- as.character(unique(dat$Condition.2))
  
  cat( "\nFitting:", condition2, "Replicate", replicate, "\n" )
  
  fit <- fit_growthcurve(dat)
   if(is.null(fit)) {
      warning(
      paste( "Growthcurver failed:", condition2, "Replicate", replicate ) )
    next
  }
  
  # Extracting parameters
  
  parameter_row <- data.frame( Condition = condition, pH = ph, Replicate = replicate,
    Condition.2 = condition2, k = fit$vals$k, n0 = fit$vals$n0, r = fit$vals$r,
    t_mid = fit$vals$t_mid, t_gen = fit$vals$t_gen,  auc_l = fit$vals$auc_l, 
    auc_e = fit$vals$auc_e,  sigma = fit$vals$sigma)
  
  all_parameters[[length(all_parameters) + 1]] <- parameter_row
  
  prediction_time <- seq( min(dat$Time.d), max(dat$Time.d), length.out = 500)
  predicted_nitrite <- fit$vals$k / (1 +((fit$vals$k - fit$vals$n0) / fit$vals$n0) 
                                     * exp( -fit$vals$r * prediction_time ) )
  prediction_df <- data.frame( Time.d = prediction_time, Predicted = predicted_nitrite)
  
  p_fit <- ggplot() + geom_point(  data = dat,  aes( x = Time.d, y = Nitrite.uM),
      colour = condition_colours[condition2], size = 3 ) +
    geom_line( data = prediction_df, aes( x = Time.d, y = Predicted ),
      colour = condition_colours[condition2], linewidth = 1.2 ) +
    labs(  title = paste0( condition2,  " | Replicate ", replicate ),
         subtitle = paste0("r = ", round(fit$vals$r, 3), "   k = ",
        round(fit$vals$k, 1), "   t_mid = ", round(fit$vals$t_mid, 2) ),
      x = "Time (days)", y = expression(Nitrite~(mu*M))) +
    plot_theme + theme(legend.position = "none")
  
  # Save
  
  file_name <- paste0(condition2, "_Rep", replicate, "_Growthcurver.png" )
  
  ggsave( filename = file.path (fit_dir,  file_name),
    plot = p_fit, width = 6, height = 5, dpi = 600 )
  
}

growth_parameters <- bind_rows( all_parameters)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# STATISTICS & PLOTS of extracted growth rates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Preparing dataset

file_path <- "C:/Users/s11jb5/Desktop/Growthcurver_Parameters.biofilm.pH.xlsx"

growth_parameters <- read_excel(
  file_path,
  sheet = "Growthcurver_Parameters.biofilm"
)

growth_parameters <- growth_parameters %>%
  mutate( Condition = factor( Condition, levels = c("free-living","cellobiose", "unamended")),
     pH = factor( pH, levels = c("5", "5.5", "6")),
    Condition.2 = factor( Condition.2, levels = c( "free-living5", "cellobiose5", "unamended5",
                                 "free-living5.5", "cellobiose5.5", "unamended5.5", 
                                 "free-living6", "cellobiose6", "unamended6" )  ) )


# Model for statistics
model_r <- lm( r ~ Condition.2 #pH * Condition
               , data = growth_parameters)
anova(model_r)

# Model Diagnostics
par(mfrow = c(2, 2))
plot(model_r)
shapiro.test(residuals(model_r)) # Shapiro-Wilk test 
leveneTest( r ~ Condition.2, data = growth_parameters) # Levene test 

# Kruskal-Wallis and Dunn test as the assumption of normality was not met
growth_parameters$Condition.2 <- factor(growth_parameters$Condition.2)

# Kruskal-Wallis test
kruskal_r <- kruskal.test( r ~ Condition.2, data = growth_parameters)
print(kruskal_r)

# Dunn's post hoc test with BH adjustment
dunn_r <- dunnTest( r ~ Condition.2, data = growth_parameters, method = "bh")
print(dunn_r)

# Extract pairwise comparison results
dunn_results <- dunn_r$res %>%
  rename( Comparison = Comparison, Z_value = Z,
    p_unadjusted = P.unadj, p_adjusted = P.adj ) %>%
  mutate(Significance = case_when( p_adjusted < 0.001 ~ "***",
      p_adjusted < 0.01  ~ "**", p_adjusted < 0.05  ~ "*", TRUE ~ "ns" ))

print(dunn_results)

# Summarising data for plotting

summary_parameters <- growth_parameters %>%
  group_by( Condition, pH, Condition.2) %>%
  summarise( n = sum(!is.na(r)), mean_r = mean(r, na.rm = TRUE ),
    se_r = sd( r, na.rm = TRUE ) / sqrt(n),
    .groups = "drop" )


# Plot theme
parameter_theme <- theme_classic(base_size = 14) +
  theme( axis.title = element_text( size = 14, face = "bold"),
    axis.text = element_text( size = 12, colour = "black" ),
    axis.line = element_line(colour = "black", linewidth = 0.7),
   legend.title = element_blank(), legend.text = element_text(size = 11))

# Colours
fill_colours2 <- c("free-living"= "grey50","cellobiose" = "#b47cab", "unamended"  = "#7cb486")
colour_colours2 <- c("free-living"= "grey20", "cellobiose" = "#72426A", "unamended"  = "#2d6c25")

# Growth rate plot
ggplot( summary_parameters, aes( y = pH, x = mean_r, fill = Condition,
    colour = Condition )) +
  geom_col( position = position_dodge(width = 0.75), width = 0.65,
    linewidth = 0.7 ) +
  geom_errorbar( aes(xmin = mean_r - se_r, xmax = mean_r + se_r ),
    position = position_dodge(width = 0.75), width = 0.15, linewidth = 0.7 ) +
  geom_point( data = growth_parameters, aes( y = pH,x = r, fill = Condition,
      colour = Condition ), position = position_jitterdodge( jitter.width = 0.08,
      jitter.height = 0, dodge.width = 0.75), shape = 21, size = 2, stroke = 0.8,
    inherit.aes = FALSE ) + scale_fill_manual( values = fill_colours2) +
  scale_colour_manual( values = colour_colours2 ) +
  labs(x = "Growth rate (r)", y = "pH" ) + parameter_theme


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calculating log fold change differnces in r (growth rate)
# 1. LogFC between pH 5/6 and pH 5.5/6 in Cellobiose and Unamended
# 2. LogFC between Cellobiose/Unamended for each pH
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

growth_parameters <- growth_parameters %>%
  mutate( Condition = factor( Condition, levels = c("free-living", "unamended", "cellobiose") ),
    pH = factor( pH, levels = c("5", "5.5", "6") ), log_r = log(r) )


# Model ESTIMATED pH MEANS WITHIN EACH CONDITION
model_logr <- lm( log_r ~ Condition * pH, data = growth_parameters)
anova(model_logr)
emm_pH <- emmeans( model_logr, ~ pH | Condition)

# Establishing contrasts:5 / 6 & 5.5 / 6

pH_contrasts <- contrast( emm_pH, method = list( "5 / 6"   = c(1, 0, -1),
                        "5.5 / 6" = c(0, 1, -1)), adjust = "holm")

pH_results <- as.data.frame( summary( pH_contrasts, infer = c(TRUE, TRUE) )) %>%
  mutate( lnFC = estimate, log2FC = estimate / log(2), Fold_change = exp(estimate),
    log2FC_lower = lower.CL / log(2), log2FC_upper = upper.CL / log(2),
    FC_lower = exp(lower.CL),  FC_upper = exp(upper.CL))

print(pH_results)

condition_comparisons <- pairs(pH_contrasts, by = "contrast", adjust = "tukey")
summary(condition_comparisons, infer = TRUE)

# PLOT log2FC between pH comparison

 ggplot(pH_results, aes( x = contrast, y = log2FC, fill = Condition,
   colour = Condition )) +  
   geom_hline(  yintercept = 0,  linetype = "dashed", colour = "grey50" ) +
  geom_errorbar(  aes( ymin = log2FC_lower, ymax = log2FC_upper ),  width = 0.15,
    linewidth = 0.8, position = position_dodge(width = 0.4) ) +
  geom_point(shape = 21, size = 4, stroke = 1, position = position_dodge(width = 0.4) ) +
  scale_fill_manual( values = fill_colours2 ) +
  scale_colour_manual(values = colour_colours2 ) +
  labs( y = expression(log[2]*" fold change in r"), x = "pH comparison" ) +
   parameter_theme

# Establishing contrast of cellobiose / unamended within each pH

condition_contrasts <- contrast( emm_condition, method = list(
    "cellobiose / unamended" = c(-1, 1) ), adjust = "holm")

condition_results <- as.data.frame( summary( condition_contrasts,
    infer = c(TRUE, TRUE) )) %>%
  mutate( lnFC = estimate, log2FC = estimate / log(2),
      Fold_change = exp(estimate), log2FC_lower = lower.CL / log(2),
    log2FC_upper = upper.CL / log(2), FC_lower = exp(lower.CL),
    FC_upper = exp(upper.CL) )

condition_results %>% select( pH, contrast, log2FC, Fold_change,
    log2FC_lower, log2FC_upper, FC_lower, FC_upper, SE, df, t.ratio, p.value )

#PLOT of log2FC cellobiose/unamended

ggplot( condition_results, aes( x = pH, y = log2FC )) +
  geom_errorbar( aes( ymin = log2FC_lower, ymax = log2FC_upper),
    width = 0.15, linewidth = 0.9, colour = "grey0" ) +
  geom_point( shape = 21, size = 4.5, stroke = 1, fill = "grey60",
    colour = "grey0" ) +
   labs( y = expression(log[2]*"(cellobiose / unamended)"), x = "pH") +
  parameter_theme
