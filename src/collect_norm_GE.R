#!/bin/R
.libPaths("/mnt/gaiagpfs/users/homedirs/smartinezarbas/R/x86_64-pc-linux-gnu-library/3.4")

## instaling/loading packages
#packs <- c("tidyr", "dplyr", "tibble",  "stringr", "reshape2", "readr", "edgeR")
#for (i in 1:length(packs))
#{
#    if(packs[i] %in% rownames(installed.packages()) == FALSE)
#        {
#	      install.packages(packs[i], repos="http://cran.irsn.fr/")
#  }
#  library(packs[i], character.only = TRUE)
#}

## Define arguments
#args <- commandArgs(TRUE)
#data.dir <- args[1]

#data.dir <- "/mnt/nfs/projects/ecosystem_biology/LAO/time_series/IMP_analysis/LAO_TS/CRISPR_analysis/GeneExpression"

## list all the mt files
#files <- paste(data.dir, list.files(data.dir), sep="/")

# load table, normalize counts and print new table with the normalized counts (per population)
#for (i in 1:length(files)) {
    # read table
#    read_tsv(files[i]) %>% 
#    column_to_rownames("gene") -> expr_tb
#    y <- DGEList(counts = expr_tb, group=group)
#    j <- calcNormFactors(y)
##    estimateCommonDisp.DGEList(j)
#    cpm(y, normalized.lib.sizes = TRUE) -> normalized_counts
#
#    normalized_counts %>% 
#    as.data.frame() %>% 
#    rownames_to_column("gene") %>%
#    gather(sample, value, -gene) %>%
#    mutate(date = as.Date(sample)) %>%
#    select(-sample) %>%
#    spread(date, value) -> mt_norm_pop
# 
#    fileName <- paste(files[i],"_normCounts.tsv", sep = "") 
#    write_tsv(mt_norm_pop, fileName)
#}


#########
## Normalization of total counts per gene (featureCounts files)
##########

library(stringr)
library(dplyr)
library(readr)
library(tidyr)
library(tibble)

args = commandArgs(trailingOnly=TRUE)
samples.input = args[1]
work.dir = args[2]
input.dir = args[3]
results.dir = args[4]

samples.input = "/mnt/nfs/projects/ecosystem_biology/LAO/time_series/IMP_analysis/LAO_TS/CRISPR_analysis/metadata/sample2date.tsv"
work.dir = "/mnt/nfs/projects/ecosystem_biology/LAO/Genomes/PopulationAnalysis"
input.dir = "Calculations/GeneLevel"
results.dir = "Calculations/GeneLevel/Analysed_totalCounts"


##combined table - FeatureCounts
type <- "mt"
#indir with files per sample
indir <- "Calculations/GeneLevel"
outdir <- "Calculations/Analysed_totalCounts"

setwd(work.dir)

dir.create(file.path(outdir))

# test if there is at least one argument: if not, return an error
#if (length(args)!=3) {
#    stop("Three mandatory arguments, RESULTSDIR, REPODIR and DATADIR", call.=FALSE)
#}

##samplelist
sample2date <- read_tsv(samples.input, col_name = c("sample", "date")) %>% select(1)

samples <- as.tibble(list.files(indir)) %>% separate(value, into = "sample", sep = ".m") %>% unique()
  
##binlist
bins <- list.files("Calculations/PopulationLevel") %>% gsub("_length.txt","",.) %>% .[!grepl("_gc.txt", .)]

#get initial table with gene length
for (i in 1:dim(samples)[1]) {
    infile.length <- paste(indir,"/",samples[1,1],".",type,".gene_len",".txt",sep="")
    comb <- read_tsv(infile.length, col_names = c("GeneID", "featureCounts"))
    ifelse(i==1, tmp <- comb, tmp <- rbind(tmp, comb))
    rm(comb)
}
    combined <- tmp %>% unique()
    rm(tmp)

#aggregate table for all files: all the gene counts of all the samples are collected now
for (i in 1:dim(samples)[1]) {
    infile <- list.files(indir, pattern = "mt.annotation.featureCounts.txt", full.names = TRUE) %>% .[!grepl(".summary", .)] %>% .[i]
    print(infile)
    cov_tab <- read_tsv(infile, col_names = TRUE, comment = "#") %>% select(1,2,7)
    colnames(cov_tab) <- c("GeneID", "contig", paste(samples[i,1],type,"featureCounts",sep="."))
    combined <- cov_tab %>% select(1,3) %>% full_join(combined, .)
    # only for the isolates datasets
    bins_genes_tmp <- cov_tab %>% select(1,2) %>% separate(contig, into = "bin", sep = "_NODE_")
    ifelse(i==1, bins_genes <- bins_genes_tmp, bins_genes <- unique(rbind(bins_genes, bins_genes_tmp)))
    rm(bins_genes_tmp)
}

##########
### normalization by cpm ...
##########

##tpm for combined table
#tpm <- function(counts, lengths) {
#    rate <- counts / lengths
#    rate / sum(rate) * 1e6
#}

#add new column with bin assocation ###
###
### if featureCounts does not have the bin added to the gene ID (case of the isolates, maybe I used workflow not updated?)

combined %>%
    left_join(., bins_genes, by = "GeneID") -> combined

### if featureCounts does have the bin added to the gene ID
#combined$bin=paste(str_split_fixed(combined$GeneID,"\\_",3)[,1],str_split_fixed(combined$GeneID,"\\_",3)[,2],sep="_")
#combined=data.frame(bin=combined[,ncol(combined)],combined[,-ncol(combined)])
#keep only bins in binlist
#combined=combined[combined$bin %in% bins,]


#check for incomplete cases
xx=combined[!(complete.cases(combined)),]
#short genes <100 (and few others) have NA in all values
#only continue with complete cases
combined <- combined[(complete.cases(combined)),] %>% filter(!is.na(bin))

# #convert to tpm
# combined %>% 
#   group_by(bin) %>%
#   do(as.data.frame(tpm(.[-c(1,2,3)],.$gene_length)))-> tpm.combined
# tpm.tab=cbind(combined[c(2,3)],as.data.frame(tpm.combined))  

#nasty for-loop instead as per group normalization does not work as expected
outdf <- combined[0, -2]
for (bin in levels(as.factor(combined$bin))){

  tpms.norm <- combined[combined$bin==bin,] %>%
        rename(gene_length = featureCounts) %>%
	#as.data.frame() %>% head() %>%
	select(-bin) %>%
	mutate_at(vars(-gene_length, -GeneID), funs(. / gene_length)) %>%
	select(-2) %>%
	reshape2::melt() %>% rename(sample = variable, rate = value) %>%
	group_by(sample) %>%
	mutate(tpms = rate/sum(rate)*1e6) %>% select(-rate) %>%
	spread(sample, tpms) %>%
	mutate(bin) 
    outdf <- rbind(outdf, tpms.norm) 
}

#check
#outdf[,-c(1,3)]%>%group_by(bin)%>%summarize_all(sum)

##########
#alternatively normalize with deseq2 per bin
##########
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library(DESeq2)
#########

# prepare matrix for

gene_size <- combined[,c(1,2)] %>% dplyr::rename(gene_length = featureCounts)
cond.vec=rep("TS",length(samples$sample))
type.vec=rep("1",length(samples$sample))

DESeq2.norm <- function(Xmat,cond,type)
{
    Xmat.col = ncol(Xmat)
    Xmat.row = nrow(Xmat)
    Xmat = round(Xmat)
    storage.mode(Xmat) <- 'integer'
    colData <- data.frame(condition = cond, type = type)
    dds <- DESeqDataSetFromMatrix(countData = Xmat, colData = colData, design = ~1)
        
    colData(dds)$condition = factor(colData(dds)$condition,levels=unique(cond))
    dds <- DESeq(dds, quiet = TRUE)
    dds <- estimateSizeFactors(dds)
	  #normalize the data
    YMat <- Xmat/rep(dds@colData@listData$sizeFactor, each = (Xmat.row))	    
return(YMat)
}

#similar loop
outdf2=combined[0,] %>% dplyr::rename(gene_length = featureCounts)
outdf3=combined[0,] %>% dplyr::rename(gene_length = featureCounts)
idx=0
for (bin in levels(as.factor(combined$bin))){
    idx=idx+1
    set1 <- combined[combined$bin==bin,] %>%
	column_to_rownames("GeneID")
    print(paste0(bin,": bin ", idx," out of ",length(levels(as.factor(combined$bin)))))
    Xmat <- as.matrix(select(set1, -bin, -featureCounts))
       
    Ymat=DESeq2.norm(Xmat,cond.vec,type.vec)
    
    set1 <- set1 %>% rownames_to_column("GeneID") %>% left_join(., gene_size, by = "GeneID")
    deseq.norm=data.frame(GeneID=set1$GeneID,bin=set1$bin,Ymat,gene_length=set1$featureCounts)
        #factor in length
    lengthnorm=apply(Ymat,2,function(x) (x/set1$gene_length)*1000)
    deseq_len.norm=data.frame(GeneID=set1$GeneID,bin=set1$bin,gene_length=set1$gene_length,lengthnorm)

    outdf2 <- rbind(outdf2, deseq.norm)
    outdf3 <- rbind(outdf3, deseq_len.norm)
    #remove subset from outdf and replace by normalized df
    #outdf2=outdf2[outdf2$bin!=bin,]
    #outdf2=rbind(deseq.norm,outdf2)
    #outdf3=outdf3[outdf3$bin!=bin,]
    #outdf3=rbind(deseq_len.norm,outdf3)
}

outdf[,-c(1,3)]%>%group_by(bin)%>%summarize_all(sum)
outdf2[,-c(1,3)]%>%group_by(bin)%>%summarize_all(sum)
outdf3[,-c(1,3)]%>%group_by(bin)%>%summarize_all(sum)

##Write out tables
outfile=paste(outdir,"/","genelevel_tpm.perbin_table_allisolates.tsv",sep="")
write.table(outdf,outfile,quote=F,sep="\t",row.names = F,col.names = T)

##Write out tables
outfile=paste(outdir,"/","genelevel_deseq2.perbin_table_allisolates.tsv",sep="")
write.table(outdf2,outfile,quote=F,sep="\t",row.names = F,col.names = T)

##Write out tables
outfile=paste(outdir,"/","genelevel_deseq2_lengthnorm.perbin_table_allisolates.tsv",sep="")
write.table(outdf3,outfile,quote=F,sep="\t",row.names = F,col.names = T)


#### add name of gene 
ann.dir <- "/mnt/nfs/projects/ecosystem_biology/LAO/Genomes/PopulationAnalysis/Calculations/GeneLevel"

files <- list.files(ann.dir, pattern = ".mt.annotation.bed.txt", full.name=TRUE)

for (i in 1:length(files)){
    ann.tb_tmp <- read_tsv(files[i], col_names = F) 
    ifelse(i==1, ann.tb <- ann.tb_tmp, ann.tb <- rbind(ann.tb, ann.tb_tmp))
}
 ann.tb %>% as.data.frame() %>%
     unique() %>%
    separate(X9, into = c("id", "inference", "locus_tag", "product"), sep = ";") %>%
    separate(id, into = c("tmp", "GeneID"), sep = "ID=") %>% select(GeneID, inference, locus_tag, product) %>%
    unique() -> ann.tb
    
## merge with the expression tables
    outdf %>% full_join(., ann.tb, by = "GeneID") -> outdf_ann
    outdf2 %>% full_join(., ann.tb, by = "GeneID") -> outdf2_ann
    outdf3 %>% full_join(., ann.tb, by = "GeneID") -> outdf3_ann
## save tables
    outfile=paste(outdir,"/","genelevel_tpm.perbin_table_allisolates_ann.tsv",sep="")
    write.table(outdf_ann,outfile,quote=F,sep="\t",row.names = F,col.names = T)

    ##Write out tables
    outfile=paste(outdir,"/","genelevel_deseq2.perbin_table_allisolates_ann.tsv",sep="")
    write.table(outdf2_ann,outfile,quote=F,sep="\t",row.names = F,col.names = T)

    ##Write out tables
    outfile=paste(outdir,"/","genelevel_deseq2_lengthnorm.perbin_table_allisolates_ann.tsv",sep="")
    write.table(outdf3_ann,outfile,quote=F,sep="\t",row.names = F,col.names = T)


