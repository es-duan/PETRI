#! /bin/bash
#SBATCH --job-name=PETRI_snakemake
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=esd4@uw.edu

#SBATCH --account=biology
#SBATCH --partition=compute-hugemem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=39
#SBATCH --mem=300G
#SBATCH --time=12:00:00

#SBATCH --export=all
#SBATCH --chdir=/mmfs1/gscratch/biology/kerrlab/ESD/PETRI

# Load conda module
module load conda

# Initialize Conda for the batch shell
eval "$(conda shell.bash hook)"

# Load environment
conda activate snakemake_host

# Run pipeline
snakemake --cores 39 --use-conda --rerun-incomplete
