# fNIRS meditation
# Behavioral data visualization
# This script requires one file: "behavioral_rawdata.xlsx"
# Programmed by Feng XIAO (updated on 2025.1.2)
############################################################################################################

### Preparation
## Load required packages for analysis
package_list <- c('car','tidyr','dplyr','readxl','ggpubr','ggplot2',
                  'cowplot','e1071','patchwork','gridExtra')
lapply(package_list, require, character.only = TRUE)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### Interval-timing speed
rd_time <- read_excel('behavioral_rawdata.xlsx', sheet = 'viso_time', na = "---")
rd_time$Condition <- factor(rd_time$Condition, levels = c('Before','After'))
p_time <- ggplot(data = rd_time, aes(x=Condition, y=Mean, fill=Meditation, color=Meditation, group=Meditation)) +
  geom_line(position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = Mean - Sd / sqrt(N), ymax = Mean + Sd / sqrt(N)), 
                width = 0.2, position = position_dodge(width = 0.3)) +
  labs(x = NULL, y = 'Time estimation difference') +
  ggtitle(NULL) +
  theme(
    axis.line = element_line(colour = "black", size = 0.4),
    axis.title = element_text(size = 7, color = "black"),
    axis.text = element_text(size = 7, color = "black"),
    panel.background = element_rect(fill = "transparent"),
    legend.position = 'none'
  ) +
  scale_color_manual(values = c("Endpoint" = "#B22222", "Mindfulness" = "#4169E1")) +
  scale_fill_manual(values = c("Endpoint" = "#B22222", "Mindfulness" = "#4169E1")) +
  scale_y_continuous(expand = c(0, 0), limits = c(-4, 0), oob = scales::squish,
                     breaks = seq(-4, 0, by = 1)) +
  plot_layout(nrow = 1) +
  plot_annotation(title = '(A) Interval-timing speed',
                  theme = theme(plot.title = element_text(size = 7, color = 'black',
                                                          face = 'bold')))
ggsave("pic_timing.pdf", plot = p_time, width = 2, height = 2)

### Delay discounting
rd_logk <- read_excel('behavioral_rawdata.xlsx', sheet = 'viso_logk', na = "---")
rd_logk$Condition <- factor(rd_logk$Condition, levels = c('Before','After'))
p_logk <- ggplot(data = rd_logk, aes(x=Condition, y=Mean, fill=Meditation, color=Meditation, group=Meditation)) +
  geom_line(position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = Mean - Sd / sqrt(N), ymax = Mean + Sd / sqrt(N)), 
                width = 0.2, position = position_dodge(width = 0.3)) +
  labs(x = NULL, y = 'Log k-value') +
  ggtitle(NULL) +
  theme(
    axis.line = element_line(colour = "black", size = 0.4),
    axis.title = element_text(size = 7, color = "black"),
    axis.text = element_text(size = 7, color = "black"),
    panel.background = element_rect(fill = "transparent"),
    legend.position = 'none'
  ) +
  scale_color_manual(values = c("Endpoint" = "#B22222", "Mindfulness" = "#4169E1")) +
  scale_fill_manual(values = c("Endpoint" = "#B22222", "Mindfulness" = "#4169E1")) +
  scale_y_continuous(expand = c(0, 0), limits = c(-8, 0), oob = scales::squish,
                     breaks = seq(-8, 0, by = 2)) +
  plot_layout(nrow = 1) +
  plot_annotation(title = '(B) Delay discounting',
                  theme = theme(plot.title = element_text(size = 7, color = 'black',
                                                          face = 'bold')))
ggsave("pic_logk.pdf", plot = p_logk, width = 2, height = 2)

### Monetary allocation task
rd_money <- read_excel('behavioral_rawdata.xlsx', sheet = 'viso_money', na = "---")
rd_money$Meditation <- factor(rd_money$Meditation, levels = c('Before','Endpoint',
                                                              'Mindfulness'))
rd_money$Condition <- factor(rd_money$Condition, levels = c('High','Medium','Low'))
p_money <- ggplot(data = rd_money, aes(x=Meditation, y=Mean, fill=Condition, color=Condition, group=Condition)) +
  geom_line(position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = Mean - Sd / sqrt(N), ymax = Mean + Sd / sqrt(N)), 
                width = 0.4, position = position_dodge(width = 0.3)) +
  labs(x = NULL, y = 'Money allocation (%)') +
  ggtitle(NULL) +
  theme(
    axis.line = element_line(colour = "black", size = 0.4),
    axis.title = element_text(size = 7, color = "black"),
    axis.text = element_text(size = 7, color = "black"),
    panel.background = element_rect(fill = "transparent"),
    legend.position = 'none'
  ) +
  scale_color_manual(values = c("High" = "#B22222", "Medium" = "#6F4F28",
                                "Low" = "#4169E1")) +
  scale_fill_manual(values = c("High" = "#B22222", "Medium" = "#6F4F28",
                               "Low" = "#4169E1")) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 1),
    oob = scales::squish,
    breaks = seq(0, 1, by = 0.25),
    labels = scales::label_percent(scale = 100) 
  ) +
  plot_layout(nrow = 1) +
  plot_annotation(title = '(C) Money allocation',
                  theme = theme(plot.title = element_text(size = 7, color = 'black',
                                                          face = 'bold')))
ggsave("pic_money.pdf", plot = p_money, width = 2, height = 2)
