#!/usr/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH -c 12
#SBATCH --mem=100G
#SBATCH --mail-user=asol@uoregon.edu 
#SBATCH --mail-type=ALL


/usr/bin/time -v htseq-count --stranded=yes /home/asol/bgmp/bioinfo/Bi623/QAA/star/alignment/3_2B_control_S3_L008/Mus_musculus.GRCm39.dna.primary_assemblyAligned.out.sam ../Mus_musculus.GRCm39.112.gtf
/usr/bin/time -v htseq-count --stranded=reverse /home/asol/bgmp/bioinfo/Bi623/QAA/star/alignment/3_2B_control_S3_L008/Mus_musculus.GRCm39.dna.primary_assemblyAligned.out.sam ../Mus_musculus.GRCm39.112.gtf


/usr/bin/time -v htseq-count --stranded=yes /home/asol/bgmp/bioinfo/Bi623/QAA/star/alignment/32_4G_both_S23_L008/Mus_musculus.GRCm39.dna.primary_assemblyAligned.out.sam ../Mus_musculus.GRCm39.112.gtf
/usr/bin/time -v htseq-count --stranded=reverse /home/asol/bgmp/bioinfo/Bi623/QAA/star/alignment/32_4G_both_S23_L008/Mus_musculus.GRCm39.dna.primary_assemblyAligned.out.sam ../Mus_musculus.GRCm39.112.gtf