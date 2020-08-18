#!/bin/R
.libPaths("/mnt/gaiagpfs/users/homedirs/smartinezarbas/R/x86_64-pc-linux-gnu-library/3.4")

## instaling/loading packages
packs <- c("tidyr", "dplyr", "tibble",  "stringr", "reshape2", "readr")
for (i in 1:length(packs))
{
#    if(packs[i] %in% rownames(installed.packages()) == FALSE)
#        {
#	      install.packages(packs[i], repos="http://cran.irsn.fr/")
#  }
  library(packs[i], character.only = TRUE)
}

## Define arguments
args <- commandArgs(TRUE)
sample2date.input <- args[1]
data.dir <- args[2]
out.dir <- args[3]
populationNames.input <- args[4]
#sample2date.input <- "/mnt/nfs/projects/ecosystem_biology/LAO/time_series/IMP_analysis/LAO_TS/CRISPR_analysis/metadata/sample2date.tsv"
#data.dir <- "/mnt/nfs/projects/ecosystem_biology/LAO/Genomes/PopulationAnalysis/Calculations/GeneLevel"
#out.dir <- "/mnt/nfs/projects/ecosystem_biology/LAO/Genomes/PopulationAnalysis/GeneExpression"

samples <- c("A01", "A02", "D32", "D36", "D49")
dates <- read_tsv(sample2date.input, col_names = c("sample", "date"))

#source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR")

##########
## collection of data from the average depht of coverage per gene
##########

library (edgeR)
# I need to separate the genes of each isolate
for (i in 1:5) {
  annotation_file <- paste(data.dir, "/", samples[i], ".mt.annotation.bed.txt", sep = "")
  expression_file <- paste(data.dir, "/", samples[i], ".mt.gene_depth_avg.txt", sep = "")

  annotation_tb <- read_tsv(annotation_file, col_names = FALSE) 
  expression_tb <-read_tsv(expression_file, col_names = c("geneID", "counts"))

  # new column with isolate, new column with gene ID as it appears in the expression table
    annotation_tb %>%
        separate(X1, into = "isolate", sep = "_NODE_", remove=FALSE) %>%
        separate(X9, into = "gene", sep = ";", remove=FALSE) %>%
	separate(gene, into = c("id","geneID"), sep = "=") %>% 
	select(-id) -> newAnnTb
    #write_tsv(newAnnTb, paste(out.dir, "/", samples[i], ".mt.annotation.bedX15.txt", sep = ""))
    
    newAnnTb %>%
        select(isolate, geneID, X9) %>% unique() %>%
        inner_join(expression_tb, ., by = "geneID") -> exprTb
    #write_tsv(exprTb, paste(out.dir, "/", samples[i], ".mt.geneExpAllIsolates.txt", sep = ""))

## nested loop to generate an expression table per isolate
    isolates <- newAnnTb %>%
            select(isolate) %>% unique()

    for (j in 1:length(isolates$isolate)){
	exprTb %>%
	    filter(isolate %in% isolates$isolate[j]) %>%
	    select(geneID, counts) -> isoExprTb
	#write_tsv(isoExprTb, paste(out.dir, "/", samples[i], ".mt.geneExp.", isolates$isolate[j], ".txt", sep = ""))
    }
}

# then normalize counts by sample and isolate
# list all the isolate files

# join all TPs per isolate
files <- paste(out.dir, list.files(out.dir), sep = "/") %>%
  .[grep(".mt.geneExp.", .)] %>%
  .[grep("Isolate_", .)]

isolates <- files %>% as.tibble() %>%
    separate(value, into = c("rm","iso"), sep = ".mt.geneExp.") %>%
    separate(iso, into = c("isolate"), sep = ".txt") %>%
    select(isolate) %>% unique()

for (i in 1:84){
  file <- paste(out.dir, list.files(out.dir), sep = "/") %>%
    .[grep(".mt.geneExp.", .)] %>%
    .[grep(paste(isolates$isolate[i], ".txt", sep = ""), .)]
  for (j in 1:length(file)) {
    samp <- file[j] %>% as.tibble() %>%
        separate(value, into = c("rm","sam"), sep = "GeneExpression/") %>%
        separate(sam, into = c("sample"), sep = ".mt.geneExp.") %>%
        select(sample) %>%
	inner_join(., dates, by = "sample") 

      read_tsv(file[j]) -> isoTb
      colnames(isoTb) <- c("gene", as.character(samp$date))

      ifelse(j==1, isoSam <- isoTb, isoSam <- full_join(isoSam, isoTb) %>% unique())
    }
    # normalization
   print(sprintf("Normalization: %s", isolates$isolate[i]))
    isoSam %>%
     column_to_rownames("gene") -> expr_tb
     group <- factor(colnames(expr_tb))
     y <- DGEList(counts = expr_tb, group=group)
#     c <- calcNormFactors(y)
#     design <- model.matrix(~ group)
#     k <- estimateCommonDisp(c)
#     estimateCommonDisp.DGEList(c)
     cpm(y, normalized.lib.sizes = TRUE) -> normalized_counts
     normalized_counts %>%
     as.data.frame() %>%
     rownames_to_column("gene") %>%
     gather(sample, value, -gene) %>%
     mutate(date = as.Date(sample)) %>%
     select(-sample) %>%
     spread(date, value) -> isoTbNrmd

    file2save <- paste(out.dir, "/", isolates$isolate[i], ".mt.geneExpNorm.tsv", sep = "")
    write_tsv(isoTbNrmd, file2save)
}

##########
## collection of data from total counts  per gene
##########



