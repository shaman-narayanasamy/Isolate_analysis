
library(stringr)
library(tidyverse)
library(ggplot2)
library(data.table)
library(RColorBrewer)

##plot COG category frequencies

###TODO filenames etc need to be adapted to repo
## global vars, change these paths 
DATADIR="/media/mh/Elements/Uni_Lux/luise_paper/20172610_drep_aeromonas_luise/"
REPODIR="/home/mh/Uni_Lux/REPOS/KO_annotation/cog"

#
indir=paste(DATADIR,"COG_anno/",sep="/")
incogdir=REPODIR

##read tables
freqs=paste(DATADIR,"tables/COG_frequency_table_combined.tsv",sep="/")
cog_freqs <- read_tsv(freqs)

##plotting
# p=ggplot(cog_freqs, aes(x=lcsb_id, y=cog_cat.split )) +
#   geom_tile(aes(fill = num_assigned), color = "white") +
#   scale_fill_gradient(low = "white", high = "steelblue") 
# 
# ##heatmaply
# library(heatmaply)
# library(reshape2)
# library(webshot)
# 

#use cog_cat.split instead of name? exclude isolates
isolate.excl <- read_tsv("misc_files/isolates_excluded.txt", col_names = "isolate")
cog.freq=cog_freqs %>% filter(!(isolate %in% isolate.excl$isolate))


cog.freq$name=paste(cog.freq$cog_cat.split,cog.freq$cat_name,sep=": ")
occ=dcast(cog.freq,lcsb_id~name,value.var = "num_assigned")
rownames(occ)=occ$lcsb_id
fff=occ[-1]
# 
# heatmaply(fff,
#           plot_method="plotly",
#           file = paste(DATADIR,"heatmaply_plot_cat_nums.html",sep="/"),
#           xlab="number of genes assigned to category",
#           fontsize_row = 7,
#           fontsize_col = 7,
#           width=10,
#           height=14
#           #file = c(paste(DATADIR,"heatmaply_plot_ratioCOGCDS.html",sep="/"),
#           #         paste(DATADIR,"heatmaply_plot_ratioCOGCDS.png",sep="/"),
#           #         paste(DATADIR,"heatmaply_plot_ratioCOGCDS.pdf",sep="/"))
#           )
# webshot(url = paste(DATADIR,"heatmaply_plot_cat_nums.html",sep="/"),
#         file = paste(DATADIR,"heatmaply_plot_cat_nums.pdf",sep="/"),
#         vwidth = 1200,
#         vheight = 1500)
# 
cog.freq$name=paste(cog.freq$cog_cat.split,cog.freq$cat_name,sep="   : ")
rat=dcast(cog.freq,lcsb_id~name,value.var = "ratio_cat_totalassigned")
rownames(rat)=occ$lcsb_id
ggg=rat[-1]
# 
# heatmaply(ggg,
#           plot_method="plotly",
#           xlab="ratio num genes assigned to cat/sum genes assigned to any cat",
#           file = paste(DATADIR,"heatmaply_plot_cat_ratio.html",sep="/"),
#           fontsize_row=7,
#           fontsize_col=7
# )
# webshot(url = paste(DATADIR,"heatmaply_plot_cat_ratio.html",sep="/"),
#         file= paste(DATADIR,"heatmaply_plot_cat_ratio.pdf",sep="/"),
#         vwidth = 1000,
#         vheight = 1500)

###
library(pheatmap)
library(viridis)

#mat <- fff  ###numbers category assignments
mat <- ggg  ###ratio category assignment / all assignments

##row annotation
annot <- cog.freq %>% select(lcsb_id,Order,Family) %>% distinct() %>% as.data.frame()

rownames(annot)=annot$lcsb_id
annot <- annot %>% select(Order)
#annot <- annot %>% select(Family)

#dat <- data.frame(values = as.numeric(mat))
dat <- as.data.frame(mat) %>% gather(cog_cat,values)
ggplot(dat, aes(values)) + geom_density(bw = "SJ")

##col annotation (COG Cats)
col_groups <- substr(colnames(mat), 1, 1)
mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(mat)

##annotation colors
mat_colors <- list(Order = colorRampPalette(brewer.pal(11,"Spectral"))(length(unique(annot$Order))))
names(mat_colors$Order) <- unique(annot$Order)



pheatmap(
  mat               = mat,
  color             = viridis(12),
  #border_color      = "lightgrey",
  border_color = NA,
  show_colnames     = T,
  show_rownames     = T,
  annotation_row = annot,
  #annotation_colors=mat_colors,
  #annotation_col = mat_col,
  cellheight = 8,
  cellwidth = 15,
  # cutree_rows = 6,
  cluster_rows = T,
  cluster_cols = F,
  # drop_levels       = TRUE,
  #fontsize          = 10
  # main              = "Default Heatmap"
  filename = "misc_files/pheatmap_cog_ratios_noborder_gaps.pdf"
)

mapcogfreqs <- pheatmap(
  mat               = mat,
  color             = viridis(12),
  #border_color      = "lightgrey",
  border_color = NA,
  show_colnames     = T,
  show_rownames     = T,
  annotation_row = annot,
  #annotation_colors=mat_colors,
  #annotation_col = mat_col,
  cellheight = 8,
  cellwidth = 15,
  # cutree_rows = 6,
  cluster_rows = T,
  cluster_cols = F,
  # drop_levels       = TRUE,
  #fontsize          = 10
  # main              = "Default Heatmap"
  # filename = paste(DATADIR,"/pheatmap_cog_ratios_noborder_gaps.pdf",sep="")
)

order1= rownames(mat)[mapcogfreqs$tree_row$order]

saveRDS(file= "misc_files/row_annot.RDS", object = annot)
saveRDS(file= "misc_files/orderheatmaps.RDS", object = order1)

