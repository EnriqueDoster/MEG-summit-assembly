#!/bin/bash
#SBATCH --partition=shas
# Number of tasks needed for use case (example):
#SBATCH --ntasks=1
#
# Processors per task:
#SBATCH --qos=normal
#SBATCH --cpus-per-task=1
#SBATCH --time=20:00:00
#SBATCH --export=ALL
#SBATCH --mail-user=enriquedoster@gmail.com
#SBATCH --mail-type=ALL

### Use the multiple program flag to specify command file job
#srun --multi-prog assembly.tasks

for FILE in assembly/R*.sh; do
echo ${FILE}
sbatch ${FILE} -o ${FILE}.stdout.txt -e ${FILE}.stderr.txt
sleep 1 # pause to be kind to the scheduler
done
