#!/bin/bash

# This script requires, MACS2, TrimGalore, Bowtie and samtools
# to create fragments corresponding to ATAC seq

SRR_IDS=(SRR13288842 SRR13288843 SRR13288844 SRR13288845 \
         SRR13288838 SRR13288839 SRR13288840 SRR13288841)

BOWTIE2_INDEX="./dm6" 
THREADS=8
MACS2_OUT="macs2_output"
mkdir -p $MACS2_OUT

for id in "${SRR_IDS[@]}"; do
  echo " Processing sample $id"

  # 1. Trim (check for trimmed files jic ) 
  echo "️ Trimming adapters with shorter cutoff..."
  trim_galore --paired --length 10 ${id}_1.fastq ${id}_2.fastq
  # 2. Align (check for BAM, bowtie maps trimmed reads to the genome)
  if [[ -f ${id}.bam ]]; then
    echo " BAM file exists for $id, skipping alignment."
  else
    echo " Aligning reads..."
    bowtie2 -x $BOWTIE2_INDEX \
      -1 ${id}_1_val_1.fq \
      -2 ${id}_2_val_2.fq \
      -X 2000 -p $THREADS --very-sensitive -k 1 | \
      samtools view -bS - | samtools sort -@ $THREADS -o ${id}.bam
    samtools index ${id}.bam
  fi

  # 3. Create fragments (check for .tsv.gz, fragments of the 
  #  genomics coordinates )
  if [[ -f ${id}_fragments.tsv.gz ]]; then
    echo " Fragment file exists for $id, skipping conversion."
  else
    echo " Creating fragment file..."
    bedtools bamtobed -bedpe -i ${id}.bam | \
      awk '$1==$4 && $6-$2 < 1000 {print $1, $2, $6, $7, $9}' OFS='\t' | \
      gzip > ${id}_fragments.tsv
  fi

  # 4. Peak calling (check for MACS2 peak output, macs finds 
  #  clusters of overlapping fragments suggesting open chroma )
  if [[ -f ${MACS2_OUT}/${id}_peaks_peaks.narrowPeak ]]; then
    echo " Peaks already called for $id, skipping MACS2."
  else
    echo " Running MACS2 for $id..."
    macs2 callpeak -t ${id}.bam \
      -f BAMPE -g dm -n ${id}_peaks \
      --nomodel --shift -100 --extsize 200 \
      --outdir $MACS2_OUT
  fi

  echo " Finished $id"
  echo ""
done

echo " All samples processed or skipped if done."
