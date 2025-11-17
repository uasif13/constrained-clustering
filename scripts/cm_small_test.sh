#!/bin/bash

#SBATCH --job-name=cm_small_shared_network_n_cluster

#SBATCH --output=/project/bader/au249/logs/%x.%j.out

#SBATCH --error=/project/bader/au249/logs/%x.%j.err

#SBATCH --partition=general

#SBATCH --qos=low

#SBATCH --account=bader

#SBATCH --nodes=1

#SBATCH --ntasks-per-node=1

#SBATCH --cpus-per-task=32

#SBATCH --time=50:00:00

#SBATCH --mem=400G



# Optional: Load any necessary modules (adjust for your cluster)

module load foss/2024a igraph/0.10.16



# Set OpenMP environment (to respect --num-processors)

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK



# Create output directories if they don’t exist

mkdir -p ./output ./cc_logs



# Run the program

srun /project/bader/au249/shared/constrained-clustering/constrained_clustering CM --edgelist ./examples/test_network.tsv  --existing-clustering ./examples/test_clustering_filtered.tsv  --algorithm leiden-cpm --resolution 0.01 --num-processors 32 --output-file ./output/output_small_shared_network_clusters --log-file ./cc_logs/output_small_shared_network_clusters.log --log-level 2

