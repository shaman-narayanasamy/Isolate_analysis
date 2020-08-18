
# packages (assuming that they are installed)
library(stringr)
library(readr)
library(tibble)
library(tidyr)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(viridis)

#TODO change to R project later
setwd("/home/mh/Uni_Lux/REPOS/Isolate_analysis")

# sample and date
read_tsv("misc_files/sample_list.txt", col_names = c("sample", "date")) %>% 
  filter(sample %in% c("A01", "A02", "D32", "D36", "D49")) -> sample2date

#old:
##on gaia: /mnt/nfs/projects/ecosystem_biology/LAO/Genomes/PopulationAnalysis/Calculations/PopulationLevel/ALL.average.pop.depth.mg.tsv
# mt.dat <- read_tsv("misc_files/ALL.average.pop.depth.mt.tsv")
# mg.dat <- read_tsv("misc_files/ALL.average.pop.depth.mg.tsv")

#new misc_files/abundance tables.Rmd
mt.dat <- read_tsv("misc_files/mg_mappedReadsByTotalByIsoLenght.tsv")
mg.dat <- read_tsv("misc_files/mt_mappedReadsByTotalByIsoLenght.tsv")

##isolate list
isolist.filt <- read_tsv("misc_files/isolates_filt.txt", col_names = c("isolate","lcsbID"))
isolate.excl <- read_tsv("misc_files/isolates_excluded.txt", col_names = "isolate")

#get order and tax annotation from cog frequency plot
order1 <- readRDS(file= "misc_files/orderheatmaps.RDS")
annot <- readRDS(file= "misc_files/row_annot.RDS")

#### filter mg 
mg.filt <- mg.dat %>%
  filter(isolate %in% isolist.filt$isolate) %>%
  filter(!(isolate %in% isolate.excl$isolate)) %>%
  left_join(isolist.filt) %>%
  select(lcsbID, `2010-10-04`, `2011-01-25`, `2011-10-05`, `2011-10-12`, `2012-01-11`) %>%
  as.data.frame()

##arange according to order
mg.filt1 <- mg.filt[match(order1, mg.filt$lcsbID),]
#remove 2nd row as Isolate_00751 is not in mg depth table
mg.filt1 <- mg.filt1[-2,]

##assign rownames
rownames(mg.filt1) = mg.filt1$lcsbID
#remove ID column
mg.filt2 = mg.filt1[,-1]

#### adjust breaks for heatmaps https://slowkow.com/notes/heatmap-tutorial/
#mat=mg.filt2
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(unlist(xs), probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
#mat_breaks <- seq(min(mat), max(mat), length.out = 10)
mat_breaks <- quantile_breaks(mg.filt2, n = 50)


### plot heatmap
mgheat <- pheatmap(
  mat               = mg.filt2,
  #color             = magma(50),
  #border_color      = "lightgrey",
  #breaks = mat_breaks,
  border_color = NA,
  show_colnames     = T,
  show_rownames     = T,
  annotation_row = annot,
  #annotation_colors=mat_colors,
  #annotation_col = mat_col,
  cellheight = 8,
  cellwidth = 30,
  # cutree_rows = 6,
  cluster_rows = F,
  cluster_cols = F,
  # drop_levels       = TRUE,
  #fontsize          = 10
  # main              = "Default Heatmap"
  filename = "misc_files/pheatmap_mgrpkm_orderleveltax_normalbreaks.pdf"
)


##### SAME FOR MT
#### filter mt 
mt.filt <- mt.dat %>%
  filter(isolate %in% isolist.filt$isolate) %>%
  filter(!(isolate %in% isolate.excl$isolate)) %>%
  left_join(isolist.filt) %>%
  select(lcsbID, `2010-10-04`, `2011-01-25`, `2011-10-05`, `2011-10-12`, `2012-01-11`) %>%
  as.data.frame()

##arange according to order
mt.filt1 <- mt.filt[match(order1, mt.filt$lcsbID),]
#remove 2nd row as Isolate_00751 is not in mg depth table
mt.filt1 <- mt.filt1[-2,]

##assign rownames
rownames(mt.filt1) = mt.filt1$lcsbID
#remove ID column
mt.filt2 = mt.filt1[,-1]

mat_breaks <- quantile_breaks(mt.filt2, n = 50)

### plot heatmap
mtheat <- pheatmap(
  mat               = mt.filt2,
  # color             = magma(50),
  # breaks= mat_breaks,
  # legend_breaks = round(mat_breaks),
  #border_color      = "lightgrey",
  border_color = NA,
  show_colnames     = T,
  show_rownames     = T,
  annotation_row = annot,
  #annotation_colors=mat_colors,
  #annotation_col = mat_col,
  cellheight = 8,
  cellwidth = 30,
  # cutree_rows = 6,
  cluster_rows = F,
  cluster_cols = F,
  # drop_levels       = TRUE,
  #fontsize          = 10
  # main              = "Default Heatmap"
  filename = "misc_files/pheatmap_mtrpkm_orderleveltax.pdf"
)




