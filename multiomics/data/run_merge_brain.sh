#!/bin/bash
#SBATCH -A <id>
#SBATCH -p node
#SBATCH -t 20:00:00
#SBATCH -n 10
#SBATCH -J cr_arc_merge_brain
#SBATCH --mail-user <email-id>
#SBATCH --mail-type=ALL
#SBATCH -e wjob-%J.err
#SBATCH -o wjob-%J.out

module load bioinfo-tools

cd ../../../multiome_data/merged

../../programs/cellranger-arc-2.0.0/cellranger-arc aggr --id=brain_aggr \
                    --reference=../cellranger_data/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
                    --csv=libraries_brain.csv \
                    --normalize=none \
                    --localcores=20 \
                    --localmem=128
