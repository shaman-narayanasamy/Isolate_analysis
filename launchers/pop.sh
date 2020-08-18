#!/bin/bash -l
#OAR -n LAO_Isolates
#OAR -l nodes=1/core=36,walltime=120

source /home/users/smartinezarbas/git/gitlab/Isolate_analysis/src/preload_modules.sh

THREADS=6 snakemake -j 6 -pfk workflow_population_analysis_isolates.done -s workflows/PopulationAnalysis
