# DESeq analysis on my local machine. 

# Set working directory
setwd("~/PAI-1_and_divergent proteases")

# label column names
# colnames = c("Input","Input","Input","Input","uPA","uPA","uPA","uPA",
#                                   "FXIIa","FXIIa","FXIIa","FXIIa","uPA","Input",
#                                   "uPA","uPA","Input","Input","uPA",
#                                  "Input",
#                                   "Input","TMPRSS2","uPA","Input","TMPRSS2","uPA","Input",
#                                   "TMPRSS2","uPA")
# 
# colnames.unique = make.unique(colnames)
# 
# 
# amp_count = c(1:12)
# 
# # Taking out Thrombin and Heparin data
# for (amp in amp_count) {
#   newname = paste("amp", amp, "_counts", sep = '')
#   newname2 = paste("amp",amp,"_counts_noXB", sep = '')
#   
#   temp = read_tsv(paste('redos/Counts_amp',amp,'.txt', sep = '')) %>% 
#     select(!'...1') %>% column_to_rownames(., var = "Variants")
#   colnames(temp) = colnames.unique
#   temp = temp %>% select(!starts_with("Thr")) %>% select(!starts_with("Hep"))
#   colnames(temp) = colnames.noHeparin.noThrombin.unique
#   assign(newname, temp) 
#   
#   temp2 = temp %>% rownames_to_column(.,var = 'ID') %>% 
#     filter(!grepl('X|B', ID)) %>% column_to_rownames(., var = "ID")
#   colnames(temp2) = colnames.noHeparin.noThrombin
#   assign(newname2, temp2)
#   
#   
#   
# }

# Import amino acid substitution counts
index = c(1:12)
all_amp_data = NULL
for (i in index) {
  temp = read_tsv(paste0("counts_trimmed copy/amp",i,"_counts.txt")) %>%  #including stop codons
    rename(Variant = rowname)
  temp2 = temp %>% filter(!grepl('X|B', Variant)) %>% column_to_rownames("Variant")
  temp = temp %>% column_to_rownames("Variant")
  name = paste0("amp",i,"_counts")
  name2 = paste0("amp",i,"_counts_noXB")
  all_amp_data = rbind(all_amp_data, temp)
  assign(name, temp)
  assign(name2, temp2)
}

colData = data.frame(condition=factor(c("Input","Input","Input","Input","uPA","uPA","uPA","uPA",
                                        "FXIIa","FXIIa","FXIIa","FXIIa","uPA","Input",
                                        "uPA","uPA","Input","Input","uPA",
                                        "Input",
                                        "Input","TMPRSS2","uPA","Input","TMPRSS2","uPA","Input",
                                        "TMPRSS2","uPA")))


countData=as.matrix(all_amp_data)
dds=DESeqDataSetFromMatrix(countData, colData, formula(~ condition))
dds$condition <- relevel(dds$condition, ref = "Input")
dds=DESeq(dds, fitType = "parametric")
condition_names = resultsNames(dds)
condition_names = condition_names[2:4]

for (condition in condition_names) {
  res = results(dds, name = condition, alpha = 0.1) %>%
    as.data.frame() %>%
    rownames_to_column() %>% rename(Variant = rowname) %>% 
    na.omit()
  write_tsv(res, paste(condition,"_total.xls", sep=""))
}

res = results(dds, alpha = 0.05)

for (condition in condition_names) {
  temp =  lfcShrink(dds, coef = condition, type = 'apeglm') %>% 
    as.data.frame() %>%
    rownames_to_column() %>% rename(Variant = rowname) %>% 
    na.omit()
  assign(paste0('apeglm_shrink_', condition), temp)
  write_tsv(temp, paste0('apeglm_shrink_', condition,'.xls'))
}  

rld <- rlog(dds, blind=FALSE)
pcaData <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

plotPCA(rld, intgroup=c("condition")) #This is the PCA of the aggregate data against all amplicons

# Ensure 'condition' is a factor with levels in the desired order
pcaData$condition <- factor(pcaData$condition, levels = c("Input", "uPA", "FXIIa", "TMPRSS2"))

ggplot(pcaData, aes(PC1, PC2, color=condition, 
                    show.legend = T)) +
  geom_point(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()+
  theme_classic()+
  xlim(-25,25)+ylim(-10,20)+
  theme(
    legend.title = element_blank(),
    axis.title.x = element_text(size = 16, face = 'bold'),
    axis.title.y = element_text(size = 16, face = 'bold'),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 12))+
  scale_color_manual(values = c("#000000", "#CC79A7", "#E69F00","#56B4E9"),
                     breaks = c("Input", "uPA", "FXIIa", "TMPRSS2")) +  # set the legend order
  geom_mark_ellipse(tol = 0.01, con.cap = 0, expand = 0.02, show.legend = F)

ggsave("AggregatePCA_pretty.pdf", device = pdf, units = "in", width = 10, height = 3)

# Perform rlog transformation
rld <- rlog(dds, blind = FALSE)

# Extract the rlog matrix
rld_mat <- assay(rld)

# Compute pairwise correlation values for all samples
rld_cor <- cor(rld_mat)

# Check the output of cor() to make sure it's correct
head(rld_cor)

# Define the color palette (for optional customization)
colors <- colorRampPalette(rev(brewer.pal(9, "BrBG")))(255)

# Plot correlation heatmap using pheatmap
pheatmap(rld_cor,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         col = colors,
         main = "Sample Correlation Heatmap Including Input")

# Look individually at each amplicon# Look individually at each amplicon
amp_count = c(1:12)

for (amp in amp_count) {
  countDataName = paste('amp',amp,'_counts_noXB', sep = '') #No stops or amber stops included
  Title = paste("Amplicon", amp, sep = ' ')
  countData = as.matrix(get(countDataName))
  colData = data.frame(condition=factor(colnames))
  
  dds=DESeqDataSetFromMatrix(countData, colData, formula(~ condition))
  dds$condition <- relevel(dds$condition, ref = "Input")
  dds=DESeq(dds)
  condition_names = resultsNames(dds)
  condition_names = condition_names[2:4]
  for (condition in condition_names) {
    res = results(dds, name = condition, alpha = 0.1) %>% as.data.frame() %>%
      rownames_to_column() %>% rename(Variant = rowname)
    write_tsv(res, paste(condition,"_amp",amp,".txt", sep=""))
  }
  rld <- rlog(dds, blind=FALSE)
  sampleDists <- dist(t(assay(rld)))
  
  vsd <- DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
  pca <- prcomp(t(assay(vsd[,])))
  row.names(pca$x) <- make.unique(row.names(pca$x))
  
  components <- pca[["x"]] #pca is the object created with prcomp()
  components <- data.frame(components)
  components$PC2 <- -components$PC2
  components$PC3 <- -components$PC3
  components = cbind(components, dds$condition) %>% rename(Enzyme = "dds$condition") 
  
  tot_explained_variance_ratio <- summary(pca)[["importance"]]['Proportion of Variance',]
  totot_explained_variance_ratio <- 100 * sum(tot_explained_variance_ratio)
  
  data3d = components %>% select(Enzyme,PC1,PC2,PC3) %>% mutate(ColorCode = as.integer(factor(Enzyme)))
  head(data3d)
  assign(paste("Amp",amp,"_data3d", sep = ""), data3d)
  
  
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(rld$condition, rld$type, sep="-")
  colnames(sampleDistMatrix) <- NULL
  
  plotPCA(rld, intgroup=c("condition"))
  
  pcaData <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  niceplot<-ggplot(pcaData, aes(PC1, PC2, color=condition, 
                                show.legend = T)) +
    geom_point(size=2) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    ggtitle(Title)+
  coord_fixed()+
    theme_classic()+
    xlim(-25,25)+ylim(-10,20)+
    theme(
      legend.title = element_blank(),
      axis.title.x = element_text(size = 16, face = 'bold'),
      axis.title.y = element_text(size = 16, face = 'bold'),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      legend.text = element_text(size = 12))+
    scale_color_manual(values = c("#000000", "#CC79A7", "#E69F00","#56B4E9"),
                       breaks = c("Input", "uPA", "FXIIa", "TMPRSS2")) +  # set the legend order
    geom_mark_ellipse(tol = 0.01, con.cap = 0, expand = 0.02, show.legend = F)
  # stat_ellipse(level = 0.8)
  print(niceplot)
  store = paste('ggplot_',amp,sep = "")
  assign(store,niceplot)
}

# Nice Looking Aggregate plot of P1 vs P2 for each amplicon
pca_combined = ggarrange(ggplot_1,ggplot_2,ggplot_3,ggplot_4,ggplot_5,ggplot_6,
                         ggplot_7,ggplot_8,ggplot_9,ggplot_10,ggplot_11,ggplot_12,
                         common.legend = TRUE, nrow = 3, ncol = 4)

ggsave("PCA_P1vP2_combined.pdf", plot = pca_combined, device = pdf)



