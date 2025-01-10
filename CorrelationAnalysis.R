# fNIRS meditation
# Correlation analysis between HbO and delay discounting change
# This script requires one file: "correlation_rawdata.xlsx"
# Programmed by Feng XIAO (updated on 2025.1.5)
############################################################################################################

### Preparation
## Load required packages for analysis
package_list <- c('car','tidyr','dplyr','readxl','effsize','e1071','lmtest','mediation',
                  'ggplot2','patchwork')
lapply(package_list, require, character.only = TRUE)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### Data input
rd_endpoint <- read_excel('correlation_rawdata.xlsx', sheet = 'CH22_endpoint', na = "---")
rd_mindfulness <- read_excel('correlation_rawdata.xlsx', sheet = 'CH22_mindfulness', na = "---")

rd_endpoint <- rd_endpoint %>% 
  rename(HbO = `HbO amplitude`, Logk = `Log k change abs`)
rd_mindfulness <- rd_mindfulness %>% 
  rename(HbO = `HbO amplitude`, Logk = `Log k change abs`)

### Endpoint-focus
cor.test(rd_endpoint$HbO, rd_endpoint$Logk, method = "pearson") #r=0.44, p=.025
plot_endpoint <- ggplot(rd_endpoint, aes(x = HbO, y = Logk)) +
  geom_point(color = "#B22222", size = 0.8, alpha = 0.8) +               
  geom_smooth(method = "lm", color = "#B22222", se = TRUE, linewidth = 1) +  
  labs(title = NULL,
       x = "HbO amplitude (uM)", y = "Log k-value change") +
  xlim(-20,20)+
  ylim(-2,2)+
  ggtitle(NULL) +
  theme(
    axis.line = element_line(colour = "black", size = 0.4),
    axis.title = element_text(size = 7, color = "black"),
    axis.text = element_text(size = 7, color = "black"),
    panel.background = element_rect(fill = "transparent"),
    legend.position = 'none'
  ) +
  plot_layout(nrow = 1) +
  plot_annotation(title = '(A) Endpoint-focus meditation: MTG-delay discounting',
                  theme = theme(plot.title = element_text(size = 7, color = 'black',
                                                          face = 'bold')))
ggsave("MTG-k_endpoint.pdf", plot = plot_endpoint, width = 2.5, height = 2.5)

### Mindfulness
cor.test(rd_mindfulness$HbO, rd_mindfulness$Logk, method = "pearson") #r=0.27, p=.177
plot_mindfulness <- ggplot(rd_mindfulness, aes(x = HbO, y = Logk)) +
  geom_point(color = "#4169E1", size = 0.8, alpha = 0.8) +               
  geom_smooth(method = "lm", color = "#4169E1", se = TRUE, linewidth = 1) +  
  labs(title = NULL,
       x = "HbO amplitude (uM)", y = "Log k-value change") +
  xlim(-20,20)+
  ylim(-2,2)+
  ggtitle(NULL) +
  theme(
    axis.line = element_line(colour = "black", size = 0.4),
    axis.title = element_text(size = 7, color = "black"),
    axis.text = element_text(size = 7, color = "black"),
    panel.background = element_rect(fill = "transparent"),
    legend.position = 'none'
  ) +
  plot_layout(nrow = 1) +
  plot_annotation(title = '(B) Mindfulness meditation: MTG-delay discounting',
                  theme = theme(plot.title = element_text(size = 7, color = 'black',
                                                          face = 'bold')))
ggsave("MTG-k_mindfulness.pdf", plot = plot_mindfulness, width = 2.5, height = 2.5)