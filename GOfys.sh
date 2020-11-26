#!/usr/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --mem=4G
#SBATCH --job-name=GOfys
#SBATCH --output=GOfys%j.log
#module load lang/Python/3.6.4-intel-2018a
#module load lang/Python/3.7.0-intel-2018b
#module load lang/Miniconda3/4.4.10
source activate GOfys

#$1 = GO quers of interest, $2 = TMM normalised expression matrix. $3Trinotate annotation output

python GOfys-viking.py $1 $2 $3
