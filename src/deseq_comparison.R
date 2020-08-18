
library(data.table)
library(tidyverse)
library(imputeTS)
library(stringr)
library(DESeq2)
library(biobroom)

###globals
RESULTSDIR="/media/mh/Elements/Uni_Lux/201708_LAO/"
# REPODIR=here::here()
REPODIR="/home/mh/Uni_Lux/Niche_ecology_LAOTS"

###read data

###combined counts table
# GENELEVEL=file.path(RESULTSDIR, "PopulationAnalysis/Calculations/GeneLevel")
# counts_files <- list.files(GENELEVEL,pattern="[AD]\\d+.mt.annotation.featureCounts.txt$", full.names = TRUE)
# 
# counts_files %>%
#   map(fread,skip=1)  %>%
#   purrr::reduce(left_join, by = c("Geneid","Chr","Start","End","Strand","Length")) -> mt_counts
# colnames(mt_counts)=gsub("Mappings/([AD]\\d+)\\.mt\\.reads\\.sorted\\.bam","\\1",colnames(mt_counts))
# 
# mt.out=paste(RESULTSDIR,"PopulationAnalysis/Calculations/Analysed","mt.allbinsisolates.rawcounts.tsv",sep="/")
# write.table(mt_counts,file=mt.out,quote=F,row.names=F)


outfile=paste(RESULTSDIR,"/PopulationAnalysis/Calculations/Analysed","/","mt.allbinsisolates.rawcounts.tsv",sep="")
counts=read.table(outfile,header=T)

#load samples table
samples=read_tsv(paste(REPODIR,"aux_datafiles/sample_list.txt",sep="/"),col_names=c("sample","date"))
#remove A01 A02
# samples=samples[-c(1,2),]

#add seasons
assignSeason <- function(charv) {
  charv=as.character(charv)
  season2=vector(length=length(charv))
  for (k in 1:length(charv)){
    date=charv[k]
    parts=strsplit(date,"-")
    year=as.numeric(parts[[1]][1])
    month=as.numeric(parts[[1]][2])
    day=as.numeric(parts[[1]][3])
    if (month<=2 || (month==3 && day <= 20)){
      season2[k]="Winter"
    } else if (month<=5 || (month==6 && day <=20)) {
      season2[k]="Spring"
    } else if (month<=8 || (month==9 && day <=21)){
      season2[k]="Summer"
    } else if (month<=11 || (month==12 && day <= 21)) {
      season2[k]="Autumn"
    } else {
      season2[k]="Winter"
    }
  }
  return(season2)
}

samples %>%
  mutate(season=assignSeason(date)) -> samples

##list of isolates
inf = "/home/mh/Uni_Lux/REPOS/Isolate_analysis/misc_files/2018_03_29_List_Isolates_and_Assembly_ID.txt"
iso.list <- read_tsv(inf, col_names = c("isolateID","LCSB_ID"))


##load gff files
gff.dir= paste(RESULTSDIR,"Databases/Annotations/Isolates/",sep="/")
gff.files <- list.files(gff.dir,full.names = T) %>% 
  set_names(nm = (basename(.) %>% tools::file_path_sans_ext()))
#gff.files= gff.files[which(names(gff.files)%in%binlist)]
gffs <- map(gff.files, read_delim, delim = "\t", col_types = "ccciicccc",col_names=F) %>%
  bind_rows(.id="isolateID") %>%
  mutate(GeneID= str_replace(X9,".*ID=([^;]+);*.*","\\1")) %>%
  mutate(product=str_replace(X9,".*product=([^;]+);*.*","\\1")) %>%
  mutate(EC=str_match(X9,"[eE][cC]_[Nn]umber=(\\d+\\.[\\d-]\\.[\\d-]+\\.[\\d-]+[^;]*);")[,2]) %>%
  select(isolateID,GeneID,product,EC)


##reduce samples for isolate analysis
iso.counts <-counts %>%
  filter(grepl("Isolate",Geneid))

#add isolateID
iso.counts <- iso.counts %>%
  as.data.frame() %>%
  mutate(isolateID = str_split_fixed(Geneid, '\\|',2)[,1])

#checks
base::setdiff(iso.list$isolateID %>% sort(), gffs$isolateID %>% unique() %>% sort() )
intersect(iso.list$isolateID %>% sort(), gffs$isolateID %>% unique() %>% sort() )

base::setdiff(iso.list$isolateID %>% sort(), iso.counts$isolateID %>% unique() %>% sort() )
intersect(iso.list$isolateID %>% sort(), iso.counts$isolateID %>% unique() %>% sort() )

## subset samples
# datesofc = c("2010-10-04", "2011-01-25", "2011-10-05", "2011-10-12", "2012-01-11")
# iso.samples <- samples %>% filter(as.character(date) %in% datesofc)
seasonsofc = c("Autumn","Winter")
iso.samples <- samples %>% filter(season %in% seasonsofc)

isoofc = "Isolate_00005"
deseq_isolate <- function(isoofc) {
  iso.counts %>%
    filter(isolateID == isoofc) %>%
    select(Geneid,one_of(iso.samples$sample)) -> iso.counts.subset
  iso.counts.subset %>%
    # select(-Start,-End,-Length) %>%
    select_if(is.numeric) -> cts
  rownames(cts)= iso.counts.subset$Geneid
  
  coldata = iso.samples
  rownames(coldata)= coldata$sample
  
  ###deseq comparison
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design= ~ season)
  dds <- DESeq(dds)
  return(dds)
}
deseq_results <- function(dds){
  res <- results(dds,contrast=c("season","Autumn","Winter"),name="season_Autumn_vs_Winter")
  
  resOrdered <- res[order(res$pvalue),]
  resfin <- resOrdered %>% tidy %>% rename(gene="GeneID") %>% left_join(gffs)
  return(resfin)
}

results = tibble(isolateID = unique(iso.counts$isolateID))
results <- results %>% 
  mutate(deseq_dds = map(isolateID, function(x) deseq_isolate(x))) %>%
  mutate(deseq_res = map(deseq_dds, function(x) deseq_results(x)))

res.all <- results %>% select(-deseq_dds) %>% unnest(deseq_res)
res.sort <- res.all %>%  arrange(p.adjusted)

outf <- "/home/mh/Uni_Lux/REPOS/Isolate_analysis/misc_files/deseq_results_winter_autumn_all_tps.tsv"
write.table(res.sort, outf, sep="\t", quote= F, col.names=T, row.names = F)

#check annotation
files.dir <- "/media/mh/Elements/Uni_Lux/luise_paper/KO_annotation/"
ko.files <- list.files(pattern= ".+besthitsKO.tsv$", path=files.dir,full.names = T) %>% 
  set_names(nm = (basename(.) %>% tools::file_path_sans_ext()))
kos <- map(ko.files, read_delim, delim = "\t",col_names=T) %>%
  bind_rows(.id="isolateID") %>%
  mutate(isolateID = gsub("\\.besthitsKO","",isolateID)) %>%
  mutate(GeneID=paste(isolateID, Gene, sep="|"))
kos <- kos %>% rename(ID = "KO") 
kos <- kos %>% select(-maxScore, -hitNumber, -Gene)

res.anno <- res.sort %>% left_join(kos) %>%
  mutate(KO= gsub("ko:","",KO))

pwy_ko <- read_tsv("/home/mh/Uni_Lux/Niche_ecology_LAOTS/aux_datafiles/ko_pwy_links.tsv", col_names = T)

res.anno <- res.anno %>% left_join(pwy_ko) %>%
  left_join(gffs)

#check coverage
counts(results$deseq_dds[[1]], normalized=T)


outf <- "/home/mh/Uni_Lux/REPOS/Isolate_analysis/misc_files/deseq_results_autumn_winter_all_tps_anno.tsv"
write.table(res.sort, outf, sep="\t", quote= F, col.names=T, row.names = F)

##Filtering out spceific functions, and plotting 
library(ggrepel)

res.anno %>%
  filter(baseMean > 5, p.adjusted < 0.05) %>%
  filter(grepl("lipid",Pwy_name, ignore.case=T) | grepl("fatty acid",Pwy_name, ignore.case=T) | grepl("wax",Pwy_name, ignore.case=T)) -> res.fatty

ec_list <- c("6.2.1.3","6.2.1.1","6.2.1.13","4.2.1.17","1.1.1.35","2.3.1.16","1.3.8.-") #betaox
ec_list <- c("2.3.1.15","2.3.1.51","2.7.1.107","3.1.3.4","2.3.1.20","3.1.1.3","3.1.1.23") #tag
res.anno %>%
  filter(baseMean > 10, p.adjusted < 0.05) %>%
  filter(EC %in% ec_list) -> res.ec

res.anno %>%
  #filter(baseMean > 10, p.adjusted < 0.05) %>%
  filter(EC == "2.3.1.20") -> res.ec

res.fatty %>%
# res.ec %>%
  select(-mapID, -Pwy_name) %>% distinct() %>%
  ggplot( aes(x=isolateID, y = estimate, color=isolateID)) +
  geom_point(size=3.5) +
  geom_label_repel( aes(label= product), alpha = 0.7) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
