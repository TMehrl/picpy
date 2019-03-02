#!/bin/bash -l

#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --license=SCRATCH
#SBATCH --job-name=rss
#SBATCH --output=rss-%j-%N.out
#SBATCH --error=rss-%j-%N.err
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=tjmehrling@lbl.gov

srun -n 1 -c 4 \
	python3 $PP_PATH/picpy.py \
	rss --zeta-range -2 0 -o 4 -p 4 DATA/raw_witness_beam_*