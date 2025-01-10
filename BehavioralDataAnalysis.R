# fNIRS meditation
# Behavioral data analysis
# This script requires one file: "behavioral_rawdata.xlsx"
# Programmed by Feng XIAO (updated on 2025.1.6)
############################################################################################################

### Preparation
## Load required packages for analysis
package_list <- c('car','tidyr','dplyr','readxl','effsize','e1071','lmtest','mediation',
                  'lmtest','mediation','ggplot2','patchwork')
lapply(package_list, require, character.only = TRUE)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### Data input
rd_pretest <- read_excel('behavioral_rawdata.xlsx', sheet = 'pretest', na = "---")
rd_posttest <- read_excel('behavioral_rawdata.xlsx', sheet = 'posttest', na = "---")
rd_k <- read_excel('behavioral_rawdata.xlsx', sheet = 'k_value', na = "---")

### Demographic data analysis
## Pretest
## Age
mean(rd_pretest$Age) #20.93
sd(rd_pretest$Age) #2.09
mean((filter(rd_pretest, Gender == 1))$Age) #male: 21.89
sd((filter(rd_pretest, Gender == 1))$Age) #male: 2.26
mean((filter(rd_pretest, Gender == 2))$Age) #female: 20.65
sd((filter(rd_pretest, Gender == 2))$Age) #female: 1.99
t.test((filter(rd_pretest, Gender == 1))$Age, (filter(rd_pretest, Gender == 2))$Age,
       paired =FALSE, alternative = c("two.sided"), var.equal=FALSE,
       conf.level=0.95) #ages did not differ within males and females
##Gender
dim(filter(rd_pretest, Gender == 1))[1] #9 males
dim(filter(rd_pretest, Gender == 2))[1] #31 females
##Handedness
dim(filter(rd_pretest, Handedness == 1))[1] #left-hand: 0
dim(filter(rd_pretest, Handedness == 2))[1] #right-hand: 40
##Meditation experience
dim(filter(rd_pretest, meditation == 1))[1] #Yes: 10
dim(filter(rd_pretest, meditation == 2))[1] #No: 30
## Posttest
## Headsize
mean(rd_posttest$`Headsize (cm)`) #56.36
sd(rd_posttest$`Headsize (cm)`) #1.68

### Group comparisons
## Excluding participants with MDD history
outliers <- which(rd_posttest$`Mental disease/Long time medicine`==1) #2 participants with MDD history
Rd_pretest <- rd_pretest[-outliers,]
Rd_posttest <- rd_posttest[-outliers,]
## Meditation involvement
inv <- c(Rd_posttest$Involvement_P, Rd_posttest$Involvement_E)
mean(inv) #overall: 76.05
sd(inv) #overall: 14.97
mean(Rd_posttest$Involvement_P) #mindfulness: 76.92
sd(Rd_posttest$Involvement_P) #mindfulness: 15.55
mean(Rd_posttest$Involvement_E) #Endpoint: 75.18
sd(Rd_posttest$Involvement_E) #Endpoint: 14.52
t.test(Rd_posttest$Involvement_P, Rd_posttest$Involvement_E,paired =FALSE,
       alternative = c("two.sided"), var.equal=FALSE,
       conf.level=0.95) #meditation involvement did not differ between groups

## Interval-timing speed
# Before vs. After (Mindfulness) vs. After (Endpoint)
df_time_diff1 <- data.frame(Rd_pretest$Time_diff1, Rd_posttest$TimeP_diff1, Rd_posttest$TimeE_diff1)
colnames(df_time_diff1) <- c('Pretest','Mindfulness','Endpoint')
df_time_diff2 <- data.frame(Rd_pretest$Time_diff2, Rd_posttest$TimeP_diff2, Rd_posttest$TimeE_diff2)
colnames(df_time_diff2) <- c('Pretest','Mindfulness','Endpoint')
df_time1 <- rbind(df_time_diff1, df_time_diff2)
df_time1 <- na.omit(df_time1) #delete the missing data
df_time2 <- pivot_longer(df_time1, cols = everything(), names_to = "group", 
                         values_to = "timing_speed")
# Skewness (less than 3 is acceptable)
skewness(df_time1$Pretest) #-0.08
skewness(df_time1$Mindfulness) #-0.16
skewness(df_time1$Endpoint) #-0.66
# Lavene's tests
leveneTest(df_time2$timing_speed, factor(df_time2$group)) #equal sds
# Paired ttests with Bonferroni corrections
t.test(df_time1$Endpoint, df_time1$Pretest,
       paired=TRUE,alternative=c("two.sided"),var.equal=TRUE,conf.level=0.95,
       p.adjust.method = "bonferroni") #Endpoint: before<after, slower, p<.001
cohen.d(df_time1$Endpoint, df_time1$Pretest) #0.51
mean(df_time1$Endpoint) #-1.97
sd(df_time1$Endpoint) #1.96
mean(df_time1$Pretest) #-1
sd(df_time1$Pretest) #1.89
t.test(df_time1$Mindfulness, df_time1$Pretest,
       paired=TRUE,alternative=c("two.sided"),var.equal=TRUE,conf.level=0.95,
       p.adjust.method = "bonferroni") #Mindfulness: before<after, slower, p=.004
cohen.d(df_time1$Mindfulness, df_time1$Pretest) #0.35
mean(df_time1$Mindfulness) #-1.70
sd(df_time1$Mindfulness) #2.09
t.test(df_time1$Endpoint, df_time1$Mindfulness,
       paired=TRUE,alternative=c("two.sided"),var.equal=TRUE,conf.level=0.95,
       p.adjust.method = "bonferroni") #Endpoint=meditation, NS

## Temporal preference
# Within-group analysis: LongTerm vs. MidTerm vs. ShortTerm
# Pretest
# Proportion
df_temp_pre1 <- data.frame(Rd_pretest$Prop_high, Rd_pretest$Prop_medium, Rd_pretest$Prop_low)
colnames(df_temp_pre1) <- c('ShortTerm', 'MidTerm', 'LongTerm')
df_temp_pre2 <- pivot_longer(df_temp_pre1, cols = everything(), names_to = "group",
                          values_to = "temporal_preference")
# Skewness (less than 3 is acceptable)
skewness(df_temp_pre1$ShortTerm) #2.77
skewness(df_temp_pre1$MidTerm) #-0.18
skewness(df_temp_pre1$LongTerm) #-0.58
# Lavene's tests
leveneTest(df_temp_pre2$temporal_preference, factor(df_temp_pre2$group)) #equal sds
# Paired ttests with Bonferroni corrections
pairwise.t.test(df_temp_pre2$temporal_preference, df_temp_pre2$group,
                p.adjust.method = "bonferroni") #ShortTerm < MidTerm = LongTerm
cohen.d(df_temp_pre1$ShortTerm, df_temp_pre1$MidTerm) #0.91, ShortTerm < MidTerm
cohen.d(df_temp_pre1$ShortTerm, df_temp_pre1$LongTerm) #1.12, ShortTerm < LongTerm

# Mindfulness
# Proportion
df_temp_p1 <- data.frame(Rd_posttest$PropP_high, Rd_posttest$PropP_medium, Rd_posttest$PropP_low)
colnames(df_temp_p1) <- c('ShortTerm', 'MidTerm', 'LongTerm')
df_temp_p2 <- pivot_longer(df_temp_p1, cols = everything(), names_to = "group",
                          values_to = "temporal_preference")
# Skewness (less than 3 is acceptable)
skewness(df_temp_p1$ShorTerm) #2.20
skewness(df_temp_p1$MidTerm) #0.49
skewness(df_temp_p1$LongTerm) #-0.26
# Lavene's tests
leveneTest(df_temp_p2$temporal_preference, factor(df_temp_p2$group)) #equal sds
# Paired ttests with Bonferroni corrections
pairwise.t.test(df_temp_p2$temporal_preference, df_temp_p2$group,
                p.adjust.method = "bonferroni") #ShortTerm < MidTerm = LongTerm
cohen.d(df_temp_p1$ShortTerm, df_temp_p1$MidTerm) #0.70, ShortTerm < MidTerm
cohen.d(df_temp_p1$ShortTerm, df_temp_p1$LongTerm) #0.91, ShortTerm < LongTerm

# Endpoint
# Proportion
df_temp_e1 <- data.frame(Rd_posttest$PropE_high, Rd_posttest$PropE_medium, Rd_posttest$PropE_low)
colnames(df_temp_e1) <- c('ShortTerm', 'MidTerm', 'LongTerm')
df_temp_e2 <- pivot_longer(df_temp_e1, cols = everything(), names_to = "group",
                          values_to = "temporal_preference")
# Skewness (less than 3 is acceptable)
skewness(df_temp_e1$ShortTerm) #1.79
skewness(df_temp_e1$MidTerm) #0.72
skewness(df_temp_e1$LongTerm) #-0.17
# Lavene's tests
leveneTest(df_temp_e2$temporal_preference, factor(df_temp_e2$group)) #unequal sds
# Paired ttests with Bonferroni corrections
pairwise.t.test(df_temp_e2$temporal_preference, df_temp_e2$group,
                p.adjust.method = "bonferroni",
                pool.sd = FALSE) #ShortTerm < MidTerm = LongTerm
cohen.d(df_temp_e1$ShortTerm, df_temp_e1$MidTerm, pooled = FALSE)#0.68, ShortTerm < MidTerm
cohen.d(df_temp_e1$ShortTerm, df_temp_e1$LongTerm, pooled = FALSE) #0.75, ShortTerm < LongTerm

# Before vs. After (Mindfulness) vs. After (Endpoint)
# Long-term temporal preference
# Proportion
t.test(df_temp_p1$LongTerm, df_temp_pre1$LongTerm,
       paired=TRUE,alternative=c("two.sided"),var.equal=FALSE,conf.level=0.95,
       p.adjust.method = "bonferroni") #NS
t.test(df_temp_p1$LongTerm, df_temp_e1$LongTerm,
       paired=TRUE,alternative=c("two.sided"),var.equal=FALSE,conf.level=0.95,
       p.adjust.method = "bonferroni") #NS
t.test(df_temp_e1$LongTerm, df_temp_pre1$LongTerm,
       paired=TRUE,alternative=c("two.sided"),var.equal=FALSE,conf.level=0.95,
       p.adjust.method = "bonferroni") #NS
# Mid-term temporal preference
# Proportion
t.test(df_temp_p1$MidTerm, df_temp_pre1$MidTerm,
       paired=TRUE,alternative=c("two.sided"),var.equal=FALSE,conf.level=0.95,
       p.adjust.method = "bonferroni") #NS
t.test(df_temp_p1$MidTerm, df_temp_e1$MidTerm,
       paired=TRUE,alternative=c("two.sided"),var.equal=FALSE,conf.level=0.95,
       p.adjust.method = "bonferroni") #NS
t.test(df_temp_e1$MidTerm, df_temp_pre1$MidTerm,
       paired=TRUE,alternative=c("two.sided"),var.equal=FALSE,conf.level=0.95,
       p.adjust.method = "bonferroni") #NS
# Short-term temporal preference
# Proportion
t.test(df_temp_p1$ShortTerm, df_temp_pre1$ShortTerm,
       paired=TRUE,alternative=c("two.sided"),var.equal=FALSE,conf.level=0.95,
       p.adjust.method = "bonferroni") #NS
t.test(df_temp_p1$ShortTerm, df_temp_e1$ShortTerm,
       paired=TRUE,alternative=c("two.sided"),var.equal=FALSE,conf.level=0.95,
       p.adjust.method = "bonferroni") #NS
t.test(df_temp_e1$ShortTerm, df_temp_pre1$ShortTerm,
       paired=TRUE,alternative=c("two.sided"),var.equal=FALSE,conf.level=0.95,
       p.adjust.method = "bonferroni") #NS
# Proportion (mean + sd)
mean(df_temp_pre1$ShortTerm) #0.20
mean(df_temp_pre1$MidTerm) #0.36
mean(df_temp_pre1$LongTerm) #0.43
sd(df_temp_pre1$ShortTerm) #0.18
sd(df_temp_pre1$MidTerm) #0.18
sd(df_temp_pre1$LongTerm) #0.22

mean(df_temp_p1$ShortTerm) #0.22
mean(df_temp_p1$MidTerm) #0.36
mean(df_temp_p1$LongTerm) #0.42
sd(df_temp_p1$ShortTerm) #0.21
sd(df_temp_p1$MidTerm) #0.22
sd(df_temp_p1$LongTerm) #0.24

mean(df_temp_e1$ShortTerm) #0.21
mean(df_temp_e1$MidTerm) #0.36
mean(df_temp_e1$LongTerm) #0.42
sd(df_temp_e1$ShortTerm) #0.20
sd(df_temp_e1$MidTerm) #0.22
sd(df_temp_e1$LongTerm) #0.28

## Discounting rates (k-values)
df_k <- data.frame(rd_k$Pretest, rd_k$Mindfulness, rd_k$Endpoint)
colnames(df_k) <- c('Pretest','Mindfulness','Endpoint')
df_logk <- log(df_k) #log k values
df_logk1 <- pivot_longer(df_logk, cols = everything(), names_to = "group",
                          values_to = "discounting rate")
# Skewness (less than 3 is acceptable)
skewness(df_logk$Pretest) #-0.62
skewness(df_logk$Mindfulness) #0.93
skewness(df_logk$Endpoint) #0.90
# Lavene's tests
leveneTest(df_logk1$`discounting rate`, factor(df_logk1$group)) #equal sds
# Paired ttests with Bonferroni corrections
t.test(df_logk$Endpoint, df_logk$Pretest,
       paired=TRUE,alternative=c("two.sided"),var.equal=TRUE,conf.level=0.95,
       p.adjust.method = "bonferroni") #p=.002, End>Pretest, more SS
cohen.d(df_logk$Endpoint, df_logk$Pretest) #0.46
mean(df_logk$Endpoint) #-5.59
sd(df_logk$Endpoint) #1.77
mean(df_logk$Pretest) #-6.37
sd(df_logk$Pretest) #1.68
t.test(df_logk$Endpoint, df_logk$Mindfulness,
       paired=TRUE,alternative=c("two.sided"),var.equal=TRUE,conf.level=0.95,
       p.adjust.method = "bonferroni") #NS
t.test(df_logk$Mindfulness, df_logk$Pretest,
       paired=TRUE,alternative=c("two.sided"),var.equal=TRUE,conf.level=0.95,
       p.adjust.method = "bonferroni") #p=.002, Present>Pretest, more SS
cohen.d(df_logk$Mindfulness, df_logk$Pretest) #0.50
mean(df_logk$Mindfulness) #-5.52
sd(df_logk$Mindfulness) #1.76

## Correlation analysis between delay discounting rate changes and monetary allocation proportion changes
rd_corr1 <- data.frame(Subj = rd_pretest$SubjectNumber,
                      end_ShortTerm = rd_posttest$PropE_high-rd_pretest$Prop_high,
                      end_MidTerm = rd_posttest$PropE_medium-rd_pretest$Prop_medium,
                      end_LongTerm = rd_posttest$PropE_low-rd_pretest$Prop_low,
                      mindfulness_ShortTerm = rd_posttest$PropP_high-rd_pretest$Prop_high,
                      mindfulness_MidTerm = rd_posttest$PropP_medium-rd_pretest$Prop_medium,
                      mindfulness_LongTerm = rd_posttest$PropP_low-rd_pretest$Prop_low
)
rd_corr1 <- rd_corr1[-c(2,5),] #exclude data from MDD participants
rd_corr2 <- data.frame(end_logk = log(rd_k$Endpoint)-log(rd_k$Pretest),
                      mindfulness_logk = log(rd_k$Mindfulness)-log(rd_k$Pretest))
rd_corr <- cbind(rd_corr1,rd_corr2)
# LongTerm~logk
cor.test(rd_corr$end_LongTerm, rd_corr$end_logk, method = "pearson") #r=-0.33, p=.042
cor.test(rd_corr$mindfulness_LongTerm, rd_corr$mindfulness_logk, method = "pearson") #NS
# MidTerm~logk
cor.test(rd_corr$end_MidTerm, rd_corr$end_logk, method = "pearson") #NS
cor.test(rd_corr$mindfulness_MidTerm, rd_corr$mindfulness_logk, method = "pearson") #NS
# ShortTerm~logk
cor.test(rd_corr$end_ShortTerm, rd_corr$end_logk, method = "pearson") #r=0.55, p<.001
cor.test(rd_corr$mindfulness_ShortTerm, rd_corr$mindfulness_logk, method = "pearson") #NS

## Plotting (correlation results)
rd_corr_long1 <- data.frame(
  x = c(rd_corr$end_LongTerm, rd_corr$mindfulness_LongTerm),
  y = c(rd_corr$end_logk, rd_corr$mindfulness_logk),
  group = rep(c("Endpoint", "Mindfulness"), each = nrow(rd_corr))
)
p1<-ggplot(rd_corr_long1, aes(x = x, y = y, color = group)) +
  geom_point(size = 0.8, alpha = 0.8) +  
  geom_smooth(method = "lm", aes(fill = group), se = TRUE, linewidth = 1) +  
  labs(
    title = NULL,
    x = "For fixed deposit (post - pre)",
    y = "Log k-value (post - pre)",
    color = "Group",
    fill = "Group"
  ) +
  scale_color_manual(values = c("Endpoint" = "#B22222", "Mindfulness" = "#4169E1")) +
  scale_fill_manual(values = c("Endpoint" = "#B22222", "Mindfulness" = "#4169E1")) +
  xlim(-0.6,0.6)+
  ylim(-4,7) +
  ggtitle(NULL) +
  theme(
    axis.line = element_line(colour = "black", size = 0.4),
    axis.title = element_text(size = 7, color = "black"),
    axis.text = element_text(size = 7, color = "black"),
    panel.background = element_rect(fill = "transparent"),
    legend.position = 'none'
  ) +
  plot_layout(nrow = 1) +
  plot_annotation(title = '(A) Correlation between long-term preference and log-k value',
                  theme = theme(plot.title = element_text(size = 7, color = 'black',
                                                          face = 'bold')))
ggsave("Long-logk.pdf", plot = p1, width = 2.5, height = 2.5)

rd_corr_long2 <- data.frame(
  x = c(rd_corr$end_ShortTerm, rd_corr$mindfulness_ShortTerm),
  y = c(rd_corr$end_logk, rd_corr$mindfulness_logk),
  group = rep(c("Endpoint", "Mindfulness"), each = nrow(rd_corr))
)
p2<-ggplot(rd_corr_long2, aes(x = x, y = y, color = group)) +
  geom_point(size = 0.8, alpha = 0.8) +  
  geom_smooth(method = "lm", aes(fill = group), se = TRUE, linewidth = 1) +  
  labs(
    title = NULL,
    x = "Spend within a week (post - pre)",
    y = "Log k-value (post - pre)",
    color = "Group",
    fill = "Group"
  ) +
  xlim(-1,1)+
  ylim(-4,7) +
  scale_color_manual(values = c("Endpoint" = "#B22222", "Mindfulness" = "#4169E1")) +
  scale_fill_manual(values = c("Endpoint" = "#B22222", "Mindfulness" = "#4169E1")) +
  ggtitle(NULL) +
  theme(
    axis.line = element_line(colour = "black", size = 0.4),
    axis.title = element_text(size = 7, color = "black"),
    axis.text = element_text(size = 7, color = "black"),
    panel.background = element_rect(fill = "transparent"),
    legend.position = 'none'
  ) +
  plot_layout(nrow = 1) +
  plot_annotation(title = '(B) Correlation between short-term preference and log-k value',
                  theme = theme(plot.title = element_text(size = 7, color = 'black',
                                                          face = 'bold')))
ggsave("Short-logk.pdf", plot = p2, width = 2.5, height = 2.5)