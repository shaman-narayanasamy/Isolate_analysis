#!/bin/bash -l

cat /mnt/nfs/projects/ecosystem_biology/LAO/Genomes/Isolate_*/Assembly_2017/Analysis/results/quast/transposed_report-renamed.tsv | head -n1 > /mnt/nfs/projects/ecosystem_biology/LAO/Genomes/ALL_Isolate_Assembly_2017_stats.txt

cat /mnt/nfs/projects/ecosystem_biology/LAO/Genomes/Isolate_*/Assembly_2017/Analysis/results/quast/transposed_report-renamed.tsv | grep -v "^Assembly" >> /mnt/nfs/projects/ecosystem_biology/LAO/Genomes/ALL_Isolate_Assembly_2017_stats.txt


