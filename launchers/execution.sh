#!/bin/bash -l

### This script takes in the necessary arguments for Snakemake and attaches them into the
### appropriate Snakemake command.

source src/preload_modules.sh

date

#  CMD="SAMPLE=${1} INPUTDIR=${2} OUTPUTDIR=${2}/NewAssembly snakemake new_assembly.done --unlock"
#  echo $CMD
#
#  SAMPLE=${1} INPUTDIR=${2} OUTPUTDIR=${2}/Assembly_2017 snakemake new_assembly.done --unlock
#  
#  CMD="SAMPLE=${1} INPUTDIR=${2} OUTPUTDIR=${2}/NewAssembly snakemake -rpf new_assembly.done --rerun-incomplete"
#  echo $CMD

  #SAMPLE=${1} INPUTDIR=${2} OUTPUTDIR=${2}/Assembly_2017 snakemake -nrpf new_assembly.done --rerun-incomplete

  # Run metaquast
  THREADS=1 SAMPLE=${1} INPUTDIR=${2} OUTPUTDIR=${2}/Assembly_2017 snakemake -rpf Analysis/analysis.done --touch
  THREADS=1 SAMPLE=${1} INPUTDIR=${2} OUTPUTDIR=${2}/Assembly_2017 snakemake -rpf Analysis/rRNA/16S.fa Analysis/results/quast/transposed_report-renamed.tsv

date
