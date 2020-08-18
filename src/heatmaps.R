#!/bin/R

##### generate tables of abundance of each isolate

# packages (assuming that they are installed)
library(stringr)
library(readr)
library(tibble)
library(tidyr)
library(ggplot2)
library(dplyr)

# sample and date
read_tsv("/Users/susana.martinez/Desktop/LAO_TS_analysis/isolates/sample2date.tsv", col_names = c("sample", "date")) %>% 
  filter(sample %in% c("A01", "A02", "D32", "D36", "D49")) -> sample2date

# average of coverage depht
#Metagenomics data
data.dir <- "/Users/susana.martinez/Desktop/LAO_TS_analysis/isolates/ContigLevel"
files <-  paste(data.dir, list.files(data.dir), sep="/") %>% .[grep(".mg.contig_depth.txt", .)] 
mg <- read.table(files[1]) 
colnames(mg) <- c("contig", files[1])

for (i in 2:length(files)) {
  tmp <- read.table(files[i])
  colnames(tmp) <- c("contig", files[i])
  mg <- full_join(mg, tmp, by="contig")
}

str_split_fixed(colnames(mg)[-1], "/", n=Inf) %>% 
  as.tibble() %>% 
  dplyr::select(samples = ncol(.)) %>% 
  separate(samples, into ="sample", sep="\\.") -> tmp_colnames

colnames(mg)[2:ncol(mg)] <- tmp_colnames$sample
rm(tmp_colnames)


# aggregating the data of the isolates to calculate the average
mg %>% 
  dplyr::select(contig) %>% 
  separate(contig, into = c("isolate", "tmp"), sep='_NODE_') %>%
  mutate(contig = paste("NODE_", tmp, sep = "")) %>% 
  select(isolate) %>% 
  bind_cols(do(., data.frame(mg[grep("Isolate", mg$contig),]))) %>% 
  dplyr::select(-contig) %>% 
  gather(sample, depth, -isolate) %>% 
  group_by(isolate, sample) %>% 
  summarize(mean_depth = mean(depth, na.rm = T)) %>% 
  spread(sample, mean_depth) -> mg.isolates.pop_sample

# convert sample to date
mg.isolates.pop_sample %>% 
  gather(sample, depth, -isolate) %>% 
  full_join(., sample2date) %>% 
  group_by(date, isolate)  -> mg.isolates.pop

## save the table
mg.isolates.pop %>% select(-sample) %>% 
  spread(date, depth) %>% 
  ungroup() -> mg.iso2save
#write_tsv(mg.iso2save, "/Users/susana.martinez/Desktop/LAO_TS_analysis/isolates/ALL.average.pop.depth.mg.tsv")


## filter isolates and change name
isolates_dict <- read_tsv("/Users/susana.martinez/Desktop/LAO_TS_analysis/isolates/isolates_filt.txt", col_names = c("isolate_n", "isolate_lcsb"))

mg.iso2save %>% 
  inner_join(., isolates_dict , by = c("isolate" = "isolate_n")) %>% 
  dplyr::select(-1) %>% 
  rename(isolate = isolate_lcsb) -> mg.iso_filt
  
### heatmap of average depth per population as population abundance
TPs <- mg.iso2save %>% gather(date, mg, -isolate) %>% ungroup() %>% select(date) %>% unique()
ISOS <- mg.iso_filt %>% select(isolate) %>% arrange(isolate)
mg.iso_filt %>% 
  gather(date, mg, -isolate) %>% 
  #mutate(date = as.Date(date)) %>% 
  ggplot(aes(x=date, y=isolate)) +
    geom_tile(aes(fill = mg), color = "white") +
    scale_fill_gradient(low = "khaki", high = "blue") +
  theme_classic() +
  ggtitle("Abundance of isolates") +
  #scale_x_date(breaks = TPs$date) + 
  #scale_y_discrete(limits = MGE_order_to_plot$pspcc) +
  scale_y_discrete(limits = ISOS$isolate) +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 11)) +
  theme(axis.text.y = element_text(size = 11), 
        panel.border=element_blank()) +
  theme(strip.background = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        strip.text = element_text(size = rel(1)),
        strip.placement = "outside") +                    
  labs(y = "Isolates") + 
  labs(x = "Date") +
  theme(legend.position="right") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) -> heatmap_iso1

plotFile2save <- "/Users/susana.martinez/Desktop/LAO_TS_analysis/isolates/heatmap_filtIso_AbAb.pdf"
pdf(plotFile2save, height = 28/2.54, width = 18/2.54, useDingbats = TRUE)
print(heatmap_iso1)
dev.off()


### heatmap of the relative abundance
mg.iso_filt %>% 
  gather(date, mg, -isolate) %>% 
  group_by(date) %>% 
  mutate(total_mg = sum(mg)) %>% 
  mutate(proportion = mg*100/total_mg) %>% 
  #mutate(test = sum(proportion)) %>% 
  ungroup() %>% select(-total_mg, -mg) %>% 
  rename(mg = proportion) %>% 
  ggplot(aes(x=date, y=isolate)) +
  geom_tile(aes(fill = mg), color = "white") +
  scale_fill_gradient(low = "khaki", high = "blue") +
  theme_classic() +
  ggtitle("Relative abundance of isolates") +
  #scale_x_date(breaks = TPs$date) + 
  #scale_y_discrete(limits = MGE_order_to_plot$pspcc) +
  scale_y_discrete(limits = ISOS$isolate) +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 11)) +
  theme(axis.text.y = element_text(size = 11), 
        panel.border=element_blank()) +
  theme(strip.background = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        strip.text = element_text(size = rel(1)),
        strip.placement = "outside") +                    
  labs(y = "Isolates") + 
  labs(x = "Date") +
  theme(legend.position="right") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) -> heatmap_iso2

plotFile2save <- "/Users/susana.martinez/Desktop/LAO_TS_analysis/isolates/heatmap_filtIso_RelAb.pdf"
pdf(plotFile2save, height = 28/2.54, width = 18/2.54, useDingbats = TRUE)
print(heatmap_iso2)
dev.off()
  


#Metatranscriptomics data
data.dir <- "/Users/susana.martinez/Desktop/LAO_TS_analysis/isolates/ContigLevel"
files <-  paste(data.dir, list.files(data.dir), sep="/") %>% .[grep(".mt.contig_depth.txt", .)] 
mt <- read.table(files[1]) 
colnames(mt) <- c("contig", files[1])

for (i in 2:length(files)) {
  tmp <- read.table(files[i])
  colnames(tmp) <- c("contig", files[i])
  mt <- full_join(mt, tmp, by="contig")
}

str_split_fixed(colnames(mt)[-1], "/", n=Inf) %>% 
  as.tibble() %>% 
  dplyr::select(samples = ncol(.)) %>% 
  separate(samples, into ="sample", sep="\\.") -> tmp_colnames

colnames(mt)[2:ncol(mt)] <- tmp_colnames$sample
rm(tmp_colnames)


# aggregating the data of the isolates to calculate the average
mt %>% 
  dplyr::select(contig) %>% 
  separate(contig, into = c("isolate", "tmp"), sep='_NODE_') %>%
  mutate(contig = paste("NODE_", tmp, sep = "")) %>% 
  select(isolate) %>% 
  bind_cols(do(., data.frame(mt[grep("Isolate", mt$contig),]))) %>% 
  dplyr::select(-contig) %>% 
  gather(sample, depth, -isolate) %>% 
  group_by(isolate, sample) %>% 
  summarize(mean_depth = mean(depth, na.rm = T)) %>% 
  spread(sample, mean_depth) -> mt.isolates.pop_sample

# convert sample to date
mt.isolates.pop_sample %>% 
  gather(sample, depth, -isolate) %>% 
  full_join(., sample2date) %>% 
  group_by(date, isolate)  -> mt.isolates.pop

## save the table
mt.isolates.pop %>% select(-sample) %>% 
  spread(date, depth) %>% 
  ungroup() -> mt.iso2save
#write_tsv(mt.iso2save, "/Users/susana.martinez/Desktop/LAO_TS_analysis/isolates/ALL.average.pop.depth.mt.tsv")


## filter isolates and change name
isolates_dict <- read_tsv("/Users/susana.martinez/Desktop/LAO_TS_analysis/isolates/isolates_filt.txt", col_names = c("isolate_n", "isolate_lcsb"))

mt.iso2save %>% 
  inner_join(., isolates_dict , by = c("isolate" = "isolate_n")) %>% 
  dplyr::select(-1) %>% 
  rename(isolate = isolate_lcsb) -> mt.iso_filt

### heatmap of average depth per population as population abundance
TPs <- mt.iso2save %>% gather(date, mt, -isolate) %>% ungroup() %>% select(date) %>% unique()

mt.iso_filt %>% 
  gather(date, mt, -isolate) %>% 
  #mutate(date = as.Date(date)) %>% 
  ggplot(aes(x=date, y=isolate)) +
  geom_tile(aes(fill = mt), color = "white") +
  scale_fill_gradient(low = "lightyellow", high = "red") +
  theme_classic() +
  ggtitle("Expression of isolates") +
  #scale_x_date(breaks = TPs$date) + 
  #scale_y_discrete(limits = mtE_order_to_plot$pspcc) +
  scale_y_discrete(limits = ISOS$isolate) +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 11)) +
  theme(axis.text.y = element_text(size = 11), 
        panel.border=element_blank()) +
  theme(strip.background = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        strip.text = element_text(size = rel(1)),
        strip.placement = "outside") +                    
  labs(y = "Isolates") + 
  labs(x = "Date") +
  theme(legend.position="right") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) -> heatmap_iso1

plotFile2save <- "/Users/susana.martinez/Desktop/LAO_TS_analysis/isolates/heatmap_filtIso_AbAb_mt.pdf"
pdf(plotFile2save, height = 28/2.54, width = 18/2.54, useDingbats = TRUE)
print(heatmap_iso1)
dev.off()


## another one MT/MG

x <- mt.iso2save[,-1]/mg.iso2save[,-1] 
x$isolate <- mt.iso2save$isolate

#rownames(x) <- mt.iso2save[,1]


x %>% 
  inner_join(., isolates_dict , by = c("isolate" = "isolate_n")) %>% 
  dplyr::select(-isolate) %>% 
  rename(isolate = isolate_lcsb) -> mtmg.iso

TPs <- mt.iso2save %>% gather(date, mt, -isolate) %>% ungroup() %>% select(date) %>% unique()

mtmg.iso %>% 
  gather(date, mt, -isolate) %>% 
  #mutate(date = as.Date(date)) %>% 
  ggplot(aes(x=date, y=isolate)) +
  geom_tile(aes(fill = mt), color = "white") +
  scale_fill_gradient(low = "whitesmoke", high = "mediumpurple1") +
  theme_classic() +
  ggtitle("MT/MG ratio") +
  #scale_x_date(breaks = TPs$date) + 
  #scale_y_discrete(limits = mtE_order_to_plot$pspcc) +
  scale_y_discrete(limits = ISOS$isolate) +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 11)) +
  theme(axis.text.y = element_text(size = 11), 
        panel.border=element_blank()) +
  theme(strip.background = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        strip.text = element_text(size = rel(1)),
        strip.placement = "outside") +                    
  labs(y = "Isolates") + 
  labs(x = "Date") +
  theme(legend.position="right") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) -> heatmap_mtmg

plotFile2save <- "/Users/susana.martinez/Desktop/LAO_TS_analysis/isolates/heatmap_filtIso_ratio_mtmg.pdf"
pdf(plotFile2save, height = 28/2.54, width = 18/2.54, useDingbats = TRUE)
print(heatmap_mtmg)
dev.off()

