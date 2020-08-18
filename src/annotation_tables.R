library(stringr)
library(tidyverse)
library(ggplot2)
library(data.table)
library(readxl)

##plot COG category frequencies

## global vars, change these paths 
DATADIR="/media/mh/Elements/Uni_Lux/luise_paper/20172610_drep_aeromonas_luise/"
REPODIR="/home/mh/Uni_Lux/REPOS/KO_annotation/cog"

#
indir=paste(DATADIR,"COG_anno/",sep="/")
incogdir=REPODIR

#read in numbers of total CDS
total_cds=read.table(paste(incogdir,"total_cds_isolates.tsv",sep="/"))
colnames(total_cds)=c("isolate","sum_cds")
#read in table for COG categories and written out names
funtab=read.table(paste(incogdir,"fun2003-2014.tab",sep="/"),sep="\t",header=T)
colnames(funtab)=c("cog_cat.split","cat_name")

##prepare data, ratios of cog category occurence

filelist=list.files(indir,pattern="COGCAT.tsv",full.names = T)

read_data <- function(z){
  iso_id=str_split_fixed(basename(z),"\\.",n=3)[1,1]
  dat <- fread(z,header=T,sep="\t")
  dat$isolate=iso_id
  return(dat)
}

datalist <- lapply(filelist, read_data)

combined.tab <- rbindlist(datalist, use.names = TRUE)
colnames(combined.tab)[c(5,6)]=c("cog_id","cog_cat")

##load table with IDs and add to convert
#only isolates with lcsb id are to be used
inf = paste(REPODIR,"2018_03_29_List_Isolates_and_Assembly_ID.txt",sep="/")
names.tab <- read_tsv(inf,col_names = c("isolate","lcsb_id"))
combined.tab <- combined.tab %>% left_join(names.tab) %>% filter(grepl("[A-Z]",lcsb_id))
combined.tab$lcsb_id=gsub("(.+)","    \\1",combined.tab$lcsb_id)

combined.tab %>%
  mutate(cog_cat.split = strsplit(as.character(cog_cat), "")) %>%
  unnest() -> combi

outf=paste(DATADIR,"gene_cog_table.tsv",sep="/")
write.table(combi,outf,sep="\t",quote=F,row.names=F,col.names=T)

##add gff annotations
#load gff files
RESULTSDIR="/media/mh/Elements/Uni_Lux/201708_LAO/"
gff.dir= paste(RESULTSDIR,"Databases/Annotations/Isolates/",sep="/")
gff.files <- list.files(gff.dir,full.names = T) %>% 
  set_names(nm = (basename(.) %>% tools::file_path_sans_ext()))
#gff.files= gff.files[which(names(gff.files)%in%binlist)]
gffs <- map(gff.files, read_delim, delim = "\t", col_types = cols(),col_names=F) %>%
  bind_rows(.id="isolateID") %>%
  mutate(GeneID= str_replace(X9,".*ID=([^;]+);*.*","\\1")) %>%
  mutate(product=str_replace(X9,".*product=([^;]+);*.*","\\1")) %>%
  select(isolateID,GeneID,product)

gffs %>% 
  mutate(Gene=str_split_fixed(GeneID,"\\|",2)[,2]) -> gffs

left_join(combi,gffs) -> combi.gff

outf=paste(DATADIR,"gene_cog_table_products.tsv",sep="/")

write.table(combi.gff,outf,sep="\t",quote=F,row.names=F,col.names=T)





