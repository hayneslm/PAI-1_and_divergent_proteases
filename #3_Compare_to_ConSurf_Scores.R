#Import conSurf scores 
conserveScore = read_tsv('data_scores.txt') %>% rename(Pos = Position)

# Compare to previously published screens on the WT background using filters from original study
WTbg_0h = read_tsv("WTbg_0h_screen.txt") %>%
  filter(padj < 0.005) %>% rename(Variant = mutation)

OldNew_compare = full_join(Input_uPA_filtered, WTbg_0h, by = "Variant") %>%
  na.omit() %>% 
  rename(NewLog2 =log2FoldChange.x, OldLog2 = log2FoldChange.y)

write_tsv(OldNew_compare,"OldNew_compare.xls")

OldNew_compare_enriched = OldNew_compare %>% filter(NewLog2 >0, OldLog2 >0)
write_tsv(OldNew_compare_enriched,"OldNew_compare_enriched.xls")

ggplot(OldNew_compare, aes(x = OldLog2, y = NewLog2))+
  geom_jitter()+
  theme_classic()+
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey")

#Compare Conservation Scores (actively working on this)
all_pos = c(1:379) %>% as.numeric() %>% as_tibble() %>% rename(Pos = value)

###uPA###
uPA_pos_score = Input_uPA_filtered %>% 
  group_by(Pos) %>% summarize(sum_log = sum(log2FoldChange)) %>% 
  mutate(Pos = as.numeric(Pos))
uPA_count_score = Input_uPA_filtered %>% group_by(Pos) %>% 
  summarize(count = n()) %>% mutate(Pos = as.numeric(Pos)) %>% 
  mutate(count = as.numeric(count))
uPA_norm_count_score = full_join(uPA_pos_score, uPA_count_score, by = "Pos") %>%
  mutate(uPA_norm_count = sum_log/count)
uPA_norm_count = uPA_norm_count_score

uPA_counts_enriched = Input_uPA_filtered %>% filter(log2FoldChange>0) %>%
  group_by(Pos) %>%
  summarize(countEC = n()) %>% mutate(Pos = as.numeric(Pos)) %>%
  mutate(countEC = as.numeric(countEC))

uPA_norm_count_score_2 = full_join(uPA_counts_enriched, uPA_count_score, by = "Pos") %>%
  mutate(uPA_norm_count_2 = countEC/count)

###TMPRSS2###
TMPRSS2_pos_score = Input_TMPRSS2_filtered %>% 
  group_by(Pos) %>% summarize(sum_log = sum(log2FoldChange)) %>% 
  mutate(Pos = as.numeric(Pos))
TMPRSS2_count_score = Input_TMPRSS2_filtered %>% group_by(Pos) %>% 
  summarize(count = n()) %>% mutate(Pos = as.numeric(Pos)) %>% 
  mutate(count = as.numeric(count))
TMPRSS2_norm_count_score = full_join(TMPRSS2_pos_score, TMPRSS2_count_score, by = "Pos") %>%
  mutate(TMPRSS2_norm_count = sum_log/count)

###FXIIa###
FXIIa_pos_score = Input_FXIIa_filtered %>% 
  group_by(Pos) %>% summarize(sum_log = sum(log2FoldChange)) %>% 
  mutate(Pos = as.numeric(Pos))
FXIIa_count_score = Input_FXIIa_filtered %>% group_by(Pos) %>% 
  summarize(count = n()) %>% mutate(Pos = as.numeric(Pos)) %>% 
  mutate(count = as.numeric(count))
FXIIa_norm_count_score = full_join(FXIIa_pos_score, FXIIa_count_score, by = "Pos") %>%
  mutate(FXIIa_norm_count = sum_log/count)

##########
conservScore  = full_join(conserveScore,uPA_norm_count_score, by = "Pos") %>% 
  select(!c(sum_log,count))
conservScore = full_join(conservScore,TMPRSS2_norm_count_score, by = "Pos") %>% 
  select(!c(sum_log,count))
conservScore = full_join(conservScore,FXIIa_norm_count_score, by = "Pos") %>% 
  select(!c(sum_log,count)) 
write_tsv(conservScore, "ConserveScore.xls")

compareTemp = conservScore %>% 
  select(!TMPRSS2_norm_count) %>%
  na.omit() %>% 
  mutate(quad = if_else(uPA_norm_count>0 & FXIIa_norm_count>0,"quad_1",
                        if_else(uPA_norm_count<=0 & FXIIa_norm_count>0, "quad_2",
                                if_else(uPA_norm_count<=0 & FXIIa_norm_count<=0, "quad_3", "quad_4"))))

compare = aov(lm(score~quad, data = compareTemp))
summary(compare)

TukeyHSD(compare)


lmCompare = lm(FXIIa_norm_count~uPA_norm_count, data = compareTemp)
summary(lmCompare) #review results of linear regression

# Filter the data for desired points to label
points_to_label <- subset(compareTemp, Pos %in% c(98, 331, 251))

compareTempPlot <- ggplot(data = compareTemp, aes(x = uPA_norm_count, y = FXIIa_norm_count, 
                                                  color = score)) +
  geom_point(size = 3) +
  scale_color_gradientn(colours = rainbow(5), name = "Conservation\nScore") +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, 
              color = "black", linewidth = 1, linetype = "dashed") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16),
    axis.title.x = element_text(size = 16, face = 'bold'),
    axis.title.y = element_text(size = 16, face = 'bold'),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(hjust = 0.5, face = "bold")
  ) +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  geom_vline(xintercept = 0, linetype = 'dotted') +
  xlab("uPA Normalized Functional Score") +
  ylab("FXIIa Normalized Functional Score") +
  geom_text(data = points_to_label, aes(label = Pos), vjust = -1, hjust = 1)

# Print the plot
print(compareTempPlot)
compareTempPlot

ggsave("uPA_FXIIa_ConSurf.pdf", device = pdf)

compareTemp %>%
  group_by(quad) %>% 
  summarize(mean_score = mean(score, na.rm = TRUE),
            stdev_score = sd(score, na.rm = TRUE)) 



compareTempTMPRSS = conservScore %>% 
  select(!FXIIa_norm_count) %>%
  na.omit() %>% 
  mutate(quad = if_else(uPA_norm_count>0 & TMPRSS2_norm_count>0,"quad_1",
                        if_else(uPA_norm_count<=0 & TMPRSS2_norm_count>0, "quad_2",
                                if_else(uPA_norm_count<=0 & TMPRSS2_norm_count<=0, "quad_3", "quad_4"))))

ggplot(data = compareTempTMPRSS, aes(x = uPA_norm_count, y = TMPRSS2_norm_count, 
                                     color = score))+
  scale_color_gradientn(colours = rainbow(5))+
  geom_point(size = 3) +
  scale_color_gradientn(colours = rainbow(5), name = "Conservation\nScore") +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, 
              color = "black", linewidth = 1, linetype = "dashed") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16),
    axis.title.x = element_text(size = 16, face = 'bold'),
    axis.title.y = element_text(size = 16, face = 'bold'),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(hjust = 0.5, face = "bold")
  ) +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  geom_vline(xintercept = 0, linetype = 'dotted') +
  xlab("uPA Normalized Functional Score") +
  ylab("TMPRSS2 Normalized Functional Score")

ggsave("uPA_TMPRSS2_ConSurf.pdf", device = pdf)

compareTempTMPRSS %>%
  group_by(quad) %>% 
  summarize(mean_score = mean(score, na.rm = TRUE),
            stdev_score = sd(score, na.rm = TRUE)) 

lmCompare = lm(TMPRSS2_norm_count~uPA_norm_count, data = compareTempTMPRSS)
summary(lmCompare) #review results of linear regression

compare = aov(lm(score~quad, data = compareTempTMPRSS))
summary(compare)

TukeyHSD(compare)


###Linear Analyses### 
# Compare Normalized functional scores to conservation scores
# uPA
uPA_conservScore2 = conservScore %>% select(uPA_norm_count, score) %>%
  na.omit()
uPA_fit2 = lm(score~uPA_norm_count, data = uPA_conservScore2)
summary(uPA_fit2)

ggplot(uPA_conservScore2, aes(x = uPA_norm_count, y = score))+
  geom_point(color = "grey50", size = 2)+
  theme_classic()+
  geom_smooth(method = "lm", color = "black", se = T, linetype = "dashed")+
  theme(
    plot.title = element_text(size = 16),
    axis.title.x = element_text(size = 16, face = 'bold'),
    axis.title.y = element_text(size = 16, face = 'bold'),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14))+
  # ggtitle("uPA Normalized Score vs ConservScore")+
  xlab("Normalized Functional Mutation Score")+
  ylab("Conservation Score")

ggsave("uPA_conservScore2.pdf", device = "pdf")


#FXIIa
FXIIa_conservScore2 = conservScore %>% select(FXIIa_norm_count, score) %>%
  na.omit()
FXIIa_fit2 = lm(score~FXIIa_norm_count, data = FXIIa_conservScore2)
summary(FXIIa_fit2)

ggplot(FXIIa_conservScore2, aes(x = FXIIa_norm_count, y = score))+
  geom_point(color = "grey50", size =2)+
  theme_classic()+
  geom_smooth(method = "lm", color = "black", se = T, linetype = "dashed")+
  theme(
    plot.title = element_text(size = 16),
    axis.title.x = element_text(size = 16, face = 'bold'),
    axis.title.y = element_text(size = 16, face = 'bold'),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14))+
  # ggtitle("FXIIa Normalized Score vs ConservScore")+
  xlab("Normalized Functional Mutation Score")+
  ylab("Conservation Score")

ggsave("FXIIa_conservScore2.pdf", device = "pdf")

#TMPRSS2
TMPRSS2_conservScore2 = conservScore %>% select(TMPRSS2_norm_count, score) %>%
  na.omit()
TMPRSS2_fit2 = lm(score~TMPRSS2_norm_count, data = TMPRSS2_conservScore2)
summary(TMPRSS2_fit2)

ggplot(TMPRSS2_conservScore2, aes(x = TMPRSS2_norm_count, y = score))+
  geom_point(color = "grey50", size = 2)+
  theme_classic()+
  geom_smooth(method = "lm", color = "black", se = T, linetype = "dashed")+
  theme(
    plot.title = element_text(size = 16),
    axis.title.x = element_text(size = 16, face = 'bold'),
    axis.title.y = element_text(size = 16, face = 'bold'),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14))+
  # ggtitle("TMPRSS2 Normalized Score vs ConservScore")+
  xlab("Normalized Functional Mutation Score")+
  ylab("Conservation Score")

ggsave("TMPRSS2_conservScore2.pdf", device = "pdf")

