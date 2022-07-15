#!/bin/bash
#SBATCH --account=ctb-sgravel
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=500
#SBATCH --time=00:10:00
#SBATCH --array=1-500
#SBATCH --output=log/%x-%A_%a.out

module load r/4.1.2

path="$HOME/projects/ctb-sgravel/luke1111/cross_coalescence/"

town_pair_file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" list_of_town_pairs.csv)

townA=${town_pair_file%.csv}

echo "Job started at: `date`"

echo "within town coalescence for $SLURM_ARRAY_TASK_ID, town $town"

while read -r townB;
do
  if [ ! -f "$path/bottleneck_overlaps/overlap_${townA}_${townB}.csv" ]; then
    Rscript $path/R/compute_bottleneck_overlap.R \
    $townA $townB \
    $path/lineage_tables/ \
    $path/bottleneck_overlaps/
  fi;
done < $path/town_pairs/$town_pair_file

echo "Job finished with exit code $? at: `date`"
