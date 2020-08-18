#!/bin/bash -l
grep "" /mnt/nfs/projects/ecosystem_biology/LAO/Genomes/*/Assembly_2017/Read_stats/read-counts/preprocessed_read_counts.txt | sed -e 's:/mnt/nfs/projects/ecosystem_biology/LAO/Genomes/::g' | sed -e 's:/Assembly_2017/Read_stats/read-counts/preprocessed_read_counts.txt::' | sed -e 's/:/\t/g' > /mnt/nfs/projects/ecosystem_biology/LAO/Genomes/Read_stats/preprocessed_read_counts.txt

grep "" /mnt/nfs/projects/ecosystem_biology/LAO/Genomes/*/Assembly_2017/Read_stats/read-counts/raw_read_counts.txt | sed -e 's:/mnt/nfs/projects/ecosystem_biology/LAO/Genomes/::g' | sed -e 's:/Assembly_2017/Read_stats/read-counts/raw_read_counts.txt::' | sed -e 's/:/\t/g' > /mnt/nfs/projects/ecosystem_biology/LAO/Genomes/Read_stats/raw_read_counts.txt

paste <(cat /mnt/nfs/projects/ecosystem_biology/LAO/Genomes/Read_stats/raw_read_counts.txt) <(cat /mnt/nfs/projects/ecosystem_biology/LAO/Genomes/Read_stats/preprocessed_read_counts.txt) | cut -f1,2,5,7 > /mnt/nfs/projects/ecosystem_biology/LAO/Genomes/Read_stats/all_read_counts.txt

sed -i '1iIsolate\tUnprocessed_pairs\tPreprocessed_pairs\tPreprocessed_singletons' /mnt/nfs/projects/ecosystem_biology/LAO/Genomes/Read_stats/all_read_counts.txt
