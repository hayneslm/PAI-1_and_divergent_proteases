# This code will be for generating "fingerprints" of PAI-1 protease specificity. Run after #2

# Load Packages and Libraries
# install.packages("tidyverse");
# install.packages("reshape2");
library(tidyverse);
library(scales);
library(reshape2)
library(ggplot2)

#Input information for generating heatmaps
aa20 =  "AVLIMFYWRKHDESTNQGCP"
aa20 = unlist(strsplit(aa20, split = ""))

#use I91L sequence
PAI1seq = "VHHPPSYVAHLASDFGVRVFQQVAQASKDRNVVFSPYGVASVLAMLQLTTGGETQQQIQAAMGFKIDDKGMAPALRHLYKELMGPWNKDELSTTDAIFVQRDLKLVQGFMPHFFRLFRSTVKQVDFSEVERARFIINDWVKTHTKGMISNLLGKGAVDQLTRLVLVNALYFNGQWKTPFPDSSTHRRLFHKSDGSTVSVPMMAQTNKFNYTEFTTPDGHYYDILELPYHGDTLSMFIAAPYEKEVPLSALTNILSAQLISHWKGNMTRLPRLLVLPKFSLETEVDLRKPLENLGMTDMFRQFQADFTSLSDQEPLHVAQALQKVKIEVNESGTVASSSTAVIVSARMAPEEIIMDRPFLFVVRHNPTGTVLFMGQVMEP"
Pos = c(1:nchar(PAI1seq))
AA = unlist(strsplit(PAI1seq, split=""))
WT_coded = cbind(Pos, AA, replicate(nchar(PAI1seq),10)) %>% as.data.frame() %>% 
  rename(value = V3) %>% mutate(value = as.numeric(value)) %>%
  rename(Mut = AA) %>% 
  mutate(padj = NA) %>% 
  as_tibble()

# Heatmap for uPA selection (using all of the data but has to have at leadt 50 base mean)
Input_uPA_filt2 = Input_uPA %>% filter(baseMean>=50)

HM_data = Input_uPA_filt2 %>% select(Pos,Mut,log2FoldChange,padj) %>% 
  rename(value = log2FoldChange)
#don't have datayet about what is missing in the library...
#but not really important for this comparative analysis
HMlist = rbind(HM_data, WT_coded) %>% mutate(Pos = as.numeric(Pos)) %>% 
  filter(Pos <= 379)
HMlist     <- HMlist[order(HMlist$Pos),]
HMlist$Pos <- as.factor(HMlist$Pos)
levels(HMlist$Pos)

ggplot(HMlist, aes(x = Pos, y = Mut)) +
  geom_tile(aes(fill=value))+
  scale_x_discrete(name = "Amino Acid Position",  breaks=c(seq(0,379,50)))+
  scale_y_discrete(name = "Amino Acid", limit = rev(aa20))+
  scale_fill_gradientn(colors = c("red", "snow2", "blue", "grey30"),
                       values = rescale(c(min(HM_data$value), 0,
                                          max(HM_data$value), 10),
                                        to = c(0,1)),
                       # values = rescale(c(-5.5, 0, 6.05, 10),to = c(0,1)),
                       na.value = "ivory2") +
  theme(panel.background = element_rect(fill = 'white'),
        legend.position = "none",
        axis.text=element_text(size=10),
        axis.title = element_text(size = 14))

ggsave("uPA_fingerprint.png", device = "png", width = 11, height = 3, units = "in")


#Heatmap for FXIIa selection
Input_FXIIa_filt2 = Input_FXIIa %>% filter(baseMean>=50)

HM_data = Input_FXIIa_filt2 %>% select(Pos,Mut,log2FoldChange,padj) %>% 
  rename(value = log2FoldChange)
#don't have dataset about what is missing in the library...
#but not really important for this comparative analysis
HMlist = rbind(HM_data, WT_coded) %>% mutate(Pos = as.numeric(Pos)) %>% 
  filter(Pos <= 379)
HMlist     <- HMlist[order(HMlist$Pos),]
HMlist$Pos <- as.factor(HMlist$Pos)
# levels(HMlist$Pos)

ggplot(HMlist, aes(x = Pos, y = Mut)) +
  geom_tile(aes(fill=value))+
  scale_x_discrete(name = "Amino Acid Position",  breaks=c(seq(0,379,50)))+
  scale_y_discrete(name = "Amino Acid", limit = rev(aa20))+
  scale_fill_gradientn(colors = c("red", "snow2", "blue", "grey30"),
                       values = rescale(c(min(HM_data$value), 0,
                                          max(HM_data$value), 10),
                                        to = c(0,1)),
                       # values = rescale(c(-5.5, 0, 6.05, 10),to = c(0,1)),
                       na.value = "ivory2") +
  theme(panel.background = element_rect(fill = 'white'),
        legend.position = "none",
        axis.text=element_text(size=10),
        axis.title = element_text(size = 14))

ggsave("FXIIa_fingerprint.png", device = "png", width = 11, height = 3, units = "in")


#Heatmap for TMPRSS2 selection
Input_TMPRSS2_filt2 = Input_TMPRSS2 %>% filter(baseMean>=50)

HM_data = Input_TMPRSS2_filt2 %>% select(Pos,Mut,log2FoldChange,padj) %>% 
  rename(value = log2FoldChange)
#don't have datayet about what is missing in the library...
#but not really important for this comparative analysis
HMlist = rbind(HM_data, WT_coded) %>% mutate(Pos = as.numeric(Pos)) %>% 
  filter(Pos <= 379)
HMlist     <- HMlist[order(HMlist$Pos),]
HMlist$Pos <- as.factor(HMlist$Pos)
levels(HMlist$Pos)

ggplot(HMlist, aes(x = Pos, y = Mut)) +
  geom_tile(aes(fill=value))+
  scale_x_discrete(name = "Amino Acid Position",  breaks=c(seq(0,379,50)))+
  scale_y_discrete(name = "Amino Acid", limit = rev(aa20))+
  scale_fill_gradientn(colors = c("red", "snow2", "blue", "grey30"),
                       values = rescale(c(min(HM_data$value), 0,
                                          max(HM_data$value), 10),
                                        to = c(0,1)),
                       # values = rescale(c(-5.5, 0, 6.05, 10),to = c(0,1)),
                       na.value = "ivory2") +
  theme(panel.background = element_rect(fill = 'white'),
        legend.position = "none",
        axis.text=element_text(size=10),
        axis.title = element_text(size = 14))

ggsave("TMPRSS2_fingerprint.png", device = "png", width = 11, height = 3, units = "in")

