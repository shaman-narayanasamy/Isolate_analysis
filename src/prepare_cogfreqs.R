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
  unnest() %>%
  group_by(isolate,lcsb_id) %>%
  filter(!is.na(cog_cat.split)) %>%
  mutate(num_totalperisolate=length(cog_cat.split)) %>%
  group_by(isolate,cog_cat.split,num_totalperisolate,lcsb_id) %>%
  summarise(num_assigned=length(cog_cat.split)) -> cog_freqs.combined

cog_freqs.combined=left_join(cog_freqs.combined,total_cds)
cog_freqs.combined=left_join(cog_freqs.combined,funtab)

#add ratio
cog_freqs.combined %>% 
  mutate(lcsb_id=gsub("^\\s+","",lcsb_id)) %>% 
  mutate(lcsb_id=gsub("\\s+$","",lcsb_id)) %>%
  group_by(isolate,lcsb_id) %>%
  mutate(ratio_cat_totalassigned=num_assigned/num_totalperisolate) -> cog.freq

##add taxonomy to cog.freq
intax=paste(REPODIR,"2018_07_18_Taxonomy.xlsx",sep="/")
taxtab <- read_excel(intax,skip = 1) %>% rename(lcsb_id=`Isolate number`,isolate=Assembly)
taxtab.red <- taxtab %>%
  select(isolate,lcsb_id,Phylum,Order,Family,Genus)

cog.freq <- cog.freq %>% left_join(taxtab.red)


#write out tables
dir.create(paste(DATADIR,"tables",sep="/"))
combined=paste(DATADIR,"tables/KO_COG_annotation_table_combined_genes.tsv",sep="/")
write.table(combined.tab, file=combined,quote=F,row.names=F,col.names=T,sep="\t")

freqs=paste(DATADIR,"tables/COG_frequency_table_combined.tsv",sep="/")
write.table(cog.freq, file=freqs,quote=F,row.names=F,col.names=T,sep="\t")
