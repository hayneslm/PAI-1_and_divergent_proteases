# Compare the results of DESeq2 analyses. Make MA plots and write out results.Make a heatmap comparing common variants in each library that pass statistical threasholds.

Input_uPA = read_tsv("apeglm_shrink_condition_uPA_vs_Input.xls") %>% 
  na.omit() %>%
  mutate(score = if_else(padj<0.1, "pass", "fail"))%>% 
  mutate(temp = Variant) %>%   
  separate(temp, into = c("AA","Rest"),sep = "(?<=[A-Z])(?=[0-9])") %>% 
  separate(Rest, into = c("Pos", "Mut"), sep = "(?<=[0-9])(?=[A-Z])") %>% 
  mutate(type = if_else(AA==Mut, "wt",
                        if_else(Mut=="X", "nonsense",
                                if_else(Mut=="B", "amber","missense"))))

ggplot() + 
  scale_x_log10()+
  geom_hline(yintercept = 0, color = "black") +
  # ggtitle('Input vs uPA')+
  xlab("BaseMean Score")+
  ylab("Log2-Fold Change")+
  scale_shape_manual(values = c(1,19))+
  geom_jitter(data = subset(Input_uPA, type == "missense"), 
              aes(baseMean, log2FoldChange, color = type, shape = score), 
              stroke = 0.5, size = 1.2)+
  geom_jitter(data = subset(Input_uPA, type == "nonsense"), 
              aes(baseMean, log2FoldChange, color = type, shape = score), 
              stroke = 0.5, size = 1.2)+
  geom_jitter(data = subset(Input_uPA, type == "wt"), 
              aes(baseMean, log2FoldChange, color = type, shape = score), 
              stroke = 0.5, size = 1.2)+
  geom_jitter(data = subset(Input_uPA, type == "amber"), 
              aes(baseMean, log2FoldChange, color = type, shape = score), 
              stroke = 0.5, size = 1.2)+
  geom_vline(xintercept = 50, color = "black", linetype = "dashed")+
  theme_classic()+
  theme(
    plot.title = element_text(size = 16),
    axis.title.x = element_text(size = 16, face = 'bold'),
    axis.title.y = element_text(size = 16, face = 'bold'),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 14), # Change the legend title font
    legend.text = element_text(size = 12))+
  guides(shape = guide_none()) + # This line will remove the legend for "score"
  labs(color = "Mutation Type") # Rename 'type' legend title to 'Mutation Type'

ggsave("uPA_MA.pdf", device = pdf)

Input_FXIIa = read_tsv("apeglm_shrink_condition_FXIIa_vs_Input.xls") %>% 
  na.omit() %>%
  mutate(score = if_else(padj<0.1, "pass", "fail"))%>% 
  mutate(temp = Variant) %>%   
  separate(temp, into = c("AA","Rest"),sep = "(?<=[A-Z])(?=[0-9])") %>% 
  separate(Rest, into = c("Pos", "Mut"), sep = "(?<=[0-9])(?=[A-Z])") %>% 
  mutate(type = if_else(AA==Mut, "wt",
                        if_else(Mut=="X", "nonsense",
                                if_else(Mut=="B", "amber","missense"))))

ggplot() + 
  scale_x_log10()+
  geom_hline(yintercept = 0, color = "black") +
  # ggtitle('Input vs FXIIa')+
  xlab("Base Mean Score")+
  ylab("Log2 Fold Change")+
  scale_shape_manual(values = c(1,19))+
  geom_jitter(data = subset(Input_FXIIa, type == "missense"), 
              aes(baseMean, log2FoldChange, color = type, shape = score), 
              stroke = 0.5, size = 1.2)+
  geom_jitter(data = subset(Input_FXIIa, type == "nonsense"), 
              aes(baseMean, log2FoldChange, color = type, shape = score), 
              stroke = 0.5, size = 1.2)+
  geom_jitter(data = subset(Input_FXIIa, type == "wt"), 
              aes(baseMean, log2FoldChange, color = type, shape = score), 
              stroke = 0.5, size = 1.2)+
  geom_jitter(data = subset(Input_FXIIa, type == "amber"), 
              aes(baseMean, log2FoldChange, color = type, shape = score), 
              stroke = 0.5, size = 1.2)+
  geom_vline(xintercept = 50, color = "black", linetype = "dashed")+
  theme_classic()+
  theme(
    plot.title = element_text(size = 16),
    axis.title.x = element_text(size = 16, face = 'bold'),
    axis.title.y = element_text(size = 16, face = 'bold'),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 14), # Change the legend title font
    legend.text = element_text(size = 12))+
  guides(shape = guide_none()) + # This line will remove the legend for "score"
  labs(color = "Mutation Type") # Rename 'type' legend title to 'Mutation Type'

ggsave("FXIIa_MA.pdf", device = pdf)


Input_TMPRSS2 = read_tsv("apeglm_shrink_condition_TMPRSS2_vs_Input.xls") %>% 
  na.omit() %>%
  mutate(score = if_else(padj<0.1, "pass", "fail"))%>% 
  mutate(temp = Variant) %>%   
  separate(temp, into = c("AA","Rest"),sep = "(?<=[A-Z])(?=[0-9])") %>% 
  separate(Rest, into = c("Pos", "Mut"), sep = "(?<=[0-9])(?=[A-Z])") %>% 
  mutate(type = if_else(AA==Mut, "wt",
                        if_else(Mut=="X", "nonsense",
                                if_else(Mut=="B", "amber","missense"))))

ggplot() + 
  scale_x_log10()+
  geom_hline(yintercept = 0, color = "black") +
  # ggtitle('Input vs TMPRSS2')+
  xlab("Base Mean Score")+
  ylab("Log2 Fold Change")+
  scale_shape_manual(values = c(1,19))+
  geom_jitter(data = subset(Input_TMPRSS2, type == "missense"), 
              aes(baseMean, log2FoldChange, color = type, shape = score), 
              stroke = 0.5, size = 1.2)+
  geom_jitter(data = subset(Input_TMPRSS2, type == "nonsense"), 
              aes(baseMean, log2FoldChange, color = type, shape = score), 
              stroke = 0.5, size = 1.2)+
  geom_jitter(data = subset(Input_TMPRSS2, type == "wt"), 
              aes(baseMean, log2FoldChange, color = type, shape = score), 
              stroke = 0.5, size = 1.2)+
  geom_jitter(data = subset(Input_TMPRSS2, type == "amber"), 
              aes(baseMean, log2FoldChange, color = type, shape = score), 
              stroke = 0.5, size = 1.2)+
  geom_vline(xintercept = 50, color = "black", linetype = "dashed")+
  theme_classic()+
  theme(
    plot.title = element_text(size = 16),
    axis.title.x = element_text(size = 16, face = 'bold'),
    axis.title.y = element_text(size = 16, face = 'bold'),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 14), # Change the legend title font
    legend.text = element_text(size = 12))+
  guides(shape = guide_none()) + # This line will remove the legend for "score"
  labs(color = "Mutation Type") # Rename 'type' legend title to 'Mutation Type'


ggsave("TMPRSS2_MA.pdf", device = pdf)

# Set basemean cuttoff to 50, and only keep those variants with a padj < 0.1
Input_uPA_filtered = Input_uPA %>% filter(score == 'pass', baseMean > 50,
                                          type == "missense")
Input_FXIIa_filtered = Input_FXIIa %>% filter(score == 'pass', baseMean > 50,
                                              type == "missense")
Input_TMPRSS2_filtered = Input_TMPRSS2 %>% filter(score == 'pass', baseMean > 50,
                                                  type == "missense")
# Write out results of DESeq2 analyses
write_tsv(Input_uPA, "Input_uPA.xls")
write_tsv(Input_FXIIa, "Input_FXIIa.xls")
write_tsv(Input_TMPRSS2, "Input_uPA.xls")

# Write our results of DESeq2 analysis that are filtered for significance threasholds
write_tsv(Input_uPA_filtered, "Input_uPA_filtered.xls")
write_tsv(Input_FXIIa_filtered, "Input_FXIIa_filtered.xls")
write_tsv(Input_TMPRSS2_filtered, "Input_uPA_filtered.xls")

# Compare uPA XIIa and TMPRSS2
FXIIa_match = full_join(Input_uPA_filtered,
                        Input_FXIIa_filtered,
                        by = "Variant") %>% na.omit() %>% 
  rename(Log2_uPA = log2FoldChange.x, padj_uPA = padj.x,
         Log2_XIIa = log2FoldChange.y, padj_XIIa = padj.y) %>% 
  select(Variant, Log2_uPA, padj_uPA, Log2_XIIa, padj_XIIa)

TMPRSS2_match = full_join(FXIIa_match, Input_TMPRSS2_filtered, by = "Variant") %>% 
  na.omit() %>% rename(Log2_TMPRSS2 = log2FoldChange, padj_TMPRSS2 = padj) %>% 
  select(Variant, Log2_uPA, padj_uPA, Log2_XIIa, padj_XIIa, Log2_TMPRSS2, padj_TMPRSS2)

heatmap_data = TMPRSS2_match %>% select(!contains("padj")) %>% 
  column_to_rownames(.,"Variant") %>% as.matrix()

pheatmap(heatmap_data, cluster_rows = T, show_rownames = F,
         color=colorRampPalette(c("red", "white", "blue"))(100)) 









