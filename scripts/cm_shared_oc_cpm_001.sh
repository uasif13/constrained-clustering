#!/bin/bash

#SBATCH --job-name=cm_shared_oc_cpm_001

#SBATCH --output=/project/bader/au249/logs/%x.%j.out

#SBATCH --error=/project/bader/au249/logs/%x.%j.err

#SBATCH --partition=bigmem

#SBATCH --qos=low

#SBATCH --account=bader

#SBATCH --nodes=1

#SBATCH --ntasks-per-node=1

#SBATCH --cpus-per-task=128

#SBATCH --time=8:00:00

#SBATCH --mem=1000G



# Optional: Load any necessary modules (adjust for your cluster)

module load foss/2024a igraph/0.10.16



# Set OpenMP environment (to respect --num-processors)

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

unset SLURM_MEM_PER_CPU
unset SLURM_MEM_PER_GPU
unset SLURM_MEM_PER_NODE

# Run the program

srun /project/bader/au249/shared/constrained-clustering/constrained_clustering CM --edgelist /project/bader/md724/Data/Open_Citations/open_citations_zero_indexed.tsv --existing-clustering /project/bader/md724/Data/Open_Citations/oc_v1_leiden_cpm_001_FIXED.tsv --algorithm leiden-cpm --resolution 0.01 --num-processors 128 --output-file /project/bader/au249/shared/constrained-clustering/output/output_shared_oc.tsv --log-file /project/bader/au249/shared/constrained-clustering/cc_logs/output_shared_oc.log --log-level 2

