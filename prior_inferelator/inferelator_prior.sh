#!/bin/bash
#SBATCH --job-name=motif_prior
#SBATCH --output=motif_prior.out
#SBATCH --error=motif_prior.err
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1

# Load environment and set paths
export PYTHONPATH=$PYTHONPATH:/scratch/cjs9968/project/dan_danr/clusterin_atacseq/inferelator-prior
export OPENBLAS_NUM_THREADS=1

# Run the command
python -m inferelator_prior.network_from_motifs \
  -b GSE163697_consensus_peaks_500pb_nochr.bed \
  -m CisBPDrosophilaALL.meme \
  -f dm6/fasta/genome.fa \
  -g Drosophila_melanogaster.BDGP6.54.114.gtf \
  -o inferelator_prior.tsv \
  --species fly \
  --no_tss \
  -c 4 \
  --scan FIMO \
  -w 10000 10000 \
  --tandem_window 100

