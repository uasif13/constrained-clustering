#!/bin/bash 

#SBATCH --job-name=cm_network_simple
#SBATCH --output=/project/bader/au249/logs/%x.%j.out
#SBATCH --error=/project/bader/au249/logs/%x.%j.err
#SBATCH --partition=general
#SBATCH --qos=low
#SBATCH --account=bader
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:1:00
#SBATCH --mem-per-cpu=4000M

srun /project/bader/au249/constrained-clustering/constrained_clustering CM --edgelist ./examples/test_network_simple_1.tsv --algorithm leiden-cpm --resolution 0.2 --num-processors 1 --output-file ./output/output_clusters_network_simple --log-file ./cc_logs/output.log --log-level 2
