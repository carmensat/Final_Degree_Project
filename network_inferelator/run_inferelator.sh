#!/bin/bash
#SBATCH --job-name=inferelator
#SBATCH --output=inferelator.out
#SBATCH --error=inferelator.err
#SBATCH --time=120:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --ntasks=1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cjs9968@nyu.edu

# Set environment variables
export PYTHONPATH=$PYTHONPATH:/scratch/cjs9968/project/dan_danr/clusterin_atacseq/inferelator_env
export OPENBLAS_NUM_THREADS=16
export MKL_NUM_THREADS=16
export NUMEXPR_NUM_THREADS=16

# Run the Inferelator script
echo "Running Inferelator script..."
python lamina_inferelator.py

