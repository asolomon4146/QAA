#!/usr/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH -c 12
#SBATCH --mem=80G
#SBATCH --mail-user=$asol@uoregon.edu 
#SBATCH --mail-type=ALL

#ran in QAA environment where star is installed
#ran in QAA/star/

/usr/bin/time -v STAR --runThreadN 12 --runMode alignReads \
  --outFilterMultimapNmax 3 \
  --outSAMunmapped Within KeepPairs \
  --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
  --readFilesCommand zcat \
  --readFilesIn ~/bgmp/bioinfo/Bi623/QAA/trims/sans_adapt_32_4G_both_S23_L008_R1P_001.fastq.gz \
     ~/bgmp/bioinfo/Bi623/QAA/trims/sans_adapt_32_4G_both_S23_L008_R2P_001.fastq.gz \
  --genomeDir ~/bgmp/bioinfo/Bi623/QAA/star/ \
  --outFileNamePrefix alignment/32_4G_both_S23_L008/Mus_musculus.GRCm39.dna.primary_assembly


#Ran this script but second run overwrote the first so need to rerun without the command below:

# /usr/bin/time -v STAR --runThreadN 12 --runMode alignReads \
#   --outFilterMultimapNmax 3 \
#   --outSAMunmapped Within KeepPairs \
#   --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
#   --readFilesCommand zcat \
#   --readFilesIn ~/bgmp/bioinfo/Bi623/QAA/trims/sans_adapt_3_2B_control_S3_L008_R1P_001.fastq.gz \
#      ~/bgmp/bioinfo/Bi623/QAA/trims/sans_adapt_3_2B_control_S3_L008_R2P_001.fastq.gz \
#   --genomeDir ~/bgmp/bioinfo/Bi623/QAA/star/ \
#   --outFileNamePrefix alignment/32_4G_both_S23_L008/Mus_musculus.GRCm39.dna.primary_assembly