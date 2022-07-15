#!/bin/bash
#SBATCH --account=ctb-sgravel
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2GB
#SBATCH --time=00:03:00
#SBATCH --array=1501-1709
#SBATCH --output=log/%x-%A_%a.out

module load r/4.1.2

path="$HOME/projects/ctb-sgravel/luke1111/cross_coalescence/"

town=$(sed -n "${SLURM_ARRAY_TASK_ID}p" list_of_towns.csv)

echo "Job started at: `date`"

echo "within town coalescence for $SLURM_ARRAY_TASK_ID, town $town"

Rscript $path/R/within_town_coalescence.R $path/data/tout_balsac.csv $town $path/lineage_tables/

echo "Job finished with exit code $? at: `date`"
