#!/bin/bash -l

## may be started using
## sbatch --dependency=afterok:<jobID_A> rss_batch.sh 

#SBATCH --partition=regular
#SBATCH --nodes=1
#SBATCH --time=06:00:00
#SBATCH --license=SCRATCH
#SBATCH --job-name=rss
#SBATCH --output=rss-%j-%N.out
#SBATCH --error=rss-%j-%N.err
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=tjmehrling@lbl.gov

Nproc=12

srun -n 1 -c $Nproc \
	python3 $PP_PATH/picpy.py rss --zeta-range -2 0 -o 4 -p $Nproc DATA/raw_witness_beam_*
