#!/bin/bash -l

### This script takes in the necessary arguments for Snakemake and attaches them into the
### appropriate Snakemake command.

source src/preload_modules.sh
date

  SAMPLE=${1} INPUTDIR=${2} OUTPUTDIR=${2}/Assembly_2017 snakemake -Frps workflows/Read_stats 

date


