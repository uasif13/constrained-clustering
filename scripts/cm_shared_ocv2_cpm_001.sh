#!/bin/bash

#SBATCH --job-name=cm_shared_ocv2_cpm_001_w_clusters

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

unset SLURM_MEM_PER_CPU
unset SLURM_MEM_PER_GPU
unset SLURM_MEM_PER_NODE

# Run the program

srun /project/bader/au249/shared/constrained-clustering/constrained_clustering CM --edgelist /project/bader/md724/Data/Open_Citations_v2/open_citations_v2_ready.tsv --existing-clustering /project/bader/md724/Data/Open_Citations_v2/oc_v2_leiden_cpm_001_FIXED.tsv --algorithm leiden-cpm --resolution 0.01 --num-processors 32 --output-file /project/bader/au249/shared/constrained-clustering/output/output_shared_ocv2_cpm_001_w_clusters.tsv --log-file /project/bader/au249/shared/constrained-clustering/cc_logs/output_shared_ocv2_cpm_001_w_clusters.log --log-level 2

