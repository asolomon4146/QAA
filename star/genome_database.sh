#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH -c 12
#SBATCH --mem=100G
#SBATCH --time=0-3
#SBATCH --mail-user=asol@uoregon.edu 
#SBATCH --mail-type=ALL


STAR --runThreadN 12 --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles ../Mus_musculus.GRCm39.dna.primary_assembly.fa --sjdbGTFfile Mus_musculus.GRCm39.112.gtf