#!/usr/bin/env python

sam_file = "/home/asol/bgmp/bioinfo/Bi623/QAA/star/alignment/3_2B_control_S3_L008/Mus_musculus.GRCm39.dna.primary_assemblyAligned.out.sam"
mapped_counts: int = 0
unmapped_counts: int = 0

with open (sam_file, "r") as fh:
    for i, line in enumerate(fh):
        line = line.strip()
        mapped: bool = False
        primary: bool = False
        if line.startswith("K00337"):
            alignment = line.split()
            flag: int = int(alignment[1])

            if((flag & 4) != 4): #checks to see if its mapped by looking at the fourth bit to see if true or false
                mapped = True
            if((flag & 256) != 256): #if the 256 bit is false then the read is not secondary.
                primary = True

            if(mapped == True and primary == True):
                mapped_counts += 1
            elif(mapped == False and primary == True):
                unmapped_counts += 1

print(f'Number of mapped counts: {mapped_counts}')
print(f'Number of unmapped counts: {unmapped_counts}')