#!/bin/bash

#SBATCH --job-name=cm_shared_hyves_cpm_001

#SBATCH --output=/project/bader/au249/logs/%x.%j.out

#SBATCH --error=/project/bader/au249/logs/%x.%j.err

#SBATCH --partition=general

#SBATCH --qos=low

#SBATCH --account=bader

#SBATCH --nodes=1

#SBATCH --ntasks-per-node=1

#SBATCH --cpus-per-task=32

#SBATCH --time=10:00:00

#SBATCH --mem=400G



# Optional: Load any necessary modules (adjust for your cluster)

module load foss/2024a igraph/0.10.16



# Set OpenMP environment (to respect --num-processors)

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

unset SLURM_MEM_PER_CPU
unset SLURM_MEM_PER_GPU
unset SLURM_MEM_PER_NODE

# Run the program

srun /project/bader/au249/shared/constrained-clustering/constrained_clustering CM --edgelist /project/bader/au249/shared/constrained-clustering/network.tsv --existing-clustering /project/bader/au249/shared/constrained-clustering/leiden_0_01.clustering --algorithm leiden-cpm --resolution 0.01 --num-processors 1 --output-file /project/bader/au249/shared/constrained-clustering/output/output_shared_hyves.tsv --log-file /project/bader/au249/shared/constrained-clustering/cc_logs/output_shared_hyves.log --log-level 2

