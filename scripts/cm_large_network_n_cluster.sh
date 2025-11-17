#!/bin/bash 

#SBATCH --job-name=cm_large_network_n_cluster_001
#SBATCH --output=/project/bader/au249/logs/%x.%j.out
#SBATCH --error=/project/bader/au249/logs/%x.%j.err
#SBATCH --partition=general
#SBATCH --qos=low
#SBATCH --account=bader
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=50:00:00
#SBATCH --mem-per-cpu=128000M

srun /project/bader/au249/shared/constrained-clustering/constrained_clustering CM --edgelist /project/bader/md724/Data/Open_Citations_v2/open_citations_v2_ready.tsv --existing-clustering /project/bader/md724/Data/Open_Citations_v2/oc_v2_leiden_cpm_001_FIXED.tsv --algorithm leiden-cpm --resolution 0.01 --num-processors 32 --output-file ./output/output_large_network_clusters_001.tsv --log-file ./cc_logs/output_large_network_clusters_001.log --log-level 2
