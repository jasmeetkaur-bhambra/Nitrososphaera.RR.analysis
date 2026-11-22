# ______________________________________________________________________________
# Urea and ammonium effect on growth at different pH 
# ______________________________________________________________________________

#Loading dataset and packages
#Dataset as df = Dataset.Urea+Ammonium_at_varying_pH

library(readxl)
library(dplyr)
library(ggplot2)
library(car)
library(emmeans)
library(ggpubr)
library(tidyr)
library(rstatix)

df <- df %>%
  mutate( pH = factor(pH, levels = sort(unique(pH))), gr.day = as.numeric(gr.day) )

# Summarising dataset by mean and se
summary_df <- df %>% group_by(pH, Condition) %>%
  summarise(n = n(), mean_gr = mean(gr.day), sd_gr = sd(gr.day), se_gr = sd_gr / sqrt(n),
    .groups = "drop")
summary_df$pH<-as.factor(summary_df$pH)

# BAR PLOT
plot_theme <- theme_classic(base_size = 14) +
  theme(axis.text = element_text(colour = "black"), axis.title = element_text(colour = "black"),
    legend.title = element_blank(), legend.position = "top", panel.grid = element_blank())

ggplot(summary_df, aes(x = pH, y = mean_gr, fill = Condition, colour = Condition)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.5, colour = "grey35",
    linewidth = 0.4) +
  geom_errorbar( aes( ymin = mean_gr - se_gr, ymax = mean_gr + se_gr ),
    position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.6) +
  scale_fill_manual( values = c( "Urea" = "#c9a66b", "Ammonium" = "#6a7d96", 
                                 "Urea + Ammonium" = "#d84d3c" ) ) +
  scale_colour_manual( values = c("Urea" = "tan4", "Ammonium" = "#2c3e50",
                                  "Urea + Ammonium" = "red3")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  plot_theme

#Statistics: TWO-WAY ANOVA

model <- aov(gr.day ~ pH * Condition, data = df)
summary(model)

# MODEL ASSUMPTIONS
shapiro.test(residuals(model))

leveneTest(gr.day ~ interaction(pH, Condition), data = df, center = median)

par(mfrow = c(2, 2))
plot(model)

# As the assumption of normality was not met
# using the Kruskal-Wallis test, followed by Dunn’s post hoc test 

# Kruskal-Wallis test

kw_pH <- df %>% group_by(Condition) %>% kruskal_test(gr.day ~ pH)
print(kw_pH)

# Dunn's post hoc test
dunn_pH <- df %>% group_by(Condition) %>%
  dunn_test( gr.day ~ pH, p.adjust.method = "BH" ) %>%
  filter( group1 == "6" | group2 == "6")
print(dunn_pH)

# Condition differences within pH

kw_condition <- df %>% group_by(pH) %>% kruskal_test(gr.day ~ Condition)
print(kw_condition)

# Dunn's post hoc test
dunn_condition <- df %>% group_by(pH) %>%
  dunn_test(gr.day ~ Condition, p.adjust.method = "BH")
print(dunn_condition)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Log fold calculations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Reorganising dataset

df_fc <- df %>%
  mutate(pH = factor(as.character(pH), levels = c("5", "5.5", "6", "6.5", "7")),
    Condition = factor(Condition, levels = c("Urea", "Ammonium", "Urea + Ammonium"))) %>%
  filter(!is.na(pH), gr.day > 0)

df_fc %>%count(pH, Condition)

# log2fold change

model_fc <- lm(log2(gr.day) ~ pH * Condition, data = df_fc)

# pH contrasts within each Condition

emm_pH <- emmeans(model_fc, ~ pH | Condition)

ph_contrasts <- list("5 / 6"   = c(1, 0, -1, 0, 0), "5.5 / 6" = c(0, 1, -1, 0, 0),
  "6.5 / 6" = c(0, 0, -1, 1, 0), "7 / 6"   = c(0, 0, -1, 0, 1))

fc_emm <- contrast( emm_pH, method = ph_contrasts)

# Within-condition significance vs pH 6
# Dunnett adjustment for comparisons with reference

within_stats <- as.data.frame(
  summary(fc_emm, infer = c(TRUE, TRUE), adjust = "dunnett" )) %>%
  rename( log2FC = estimate, p.value = p.value)

print(within_stats)

# Between-condition statistics for each pH contrast

between_stats <- pairs( fc_emm, by = "contrast", adjust = "tukey")

between_results <- as.data.frame(summary( between_stats, infer = c(TRUE, TRUE)))

print(between_results)

# Plotting data

# Unadjusted 95% confidence intervals for effect-size visualization

plot_df <- as.data.frame(confint(fc_emm, adjust = "none")) %>%
  rename( log2FC = estimate) %>%
  mutate(contrast = factor( contrast, levels = c( "5 / 6", "5.5 / 6",
        "6.5 / 6", "7 / 6" ) ) )

# Log2 fold-change plot

 ggplot( plot_df, aes(  y = log2FC, x = contrast, fill = Condition,
    colour = Condition )) +
  geom_hline( yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar( aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.15, linewidth = 0.9, position = position_dodge(width = 0.6)) +
  geom_point(shape = 21, size = 4, stroke = 0.8, position = position_dodge(width = 0.6)) +
   scale_fill_manual( values = c( "Urea" = "#c9a66b", "Ammonium" = "#6a7d96", 
       "Urea + Ammonium" = "#d84d3c" ) ) +
  scale_colour_manual( values = c("Urea" = "tan4", "Ammonium" = "#2c3e50",
      "Urea + Ammonium" = "red3")) +
  labs(y = expression(log[2]~"fold change in growth rate"), x = "pH contrast" ) +
   plot_theme



