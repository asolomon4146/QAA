#!/usr/bin/env python

import bioinfo
import argparse
import matplotlib.pyplot as plt
import gzip
import numpy as np

def get_args():
    parser = argparse.ArgumentParser(description="A program to demultiplex")
    parser.add_argument("-b", "--barcode_list_file", help="Your file containing the list of barcodes is:", required=True, type = str)
    parser.add_argument("-r1", "--read1", help="Your read 1 file of type fastq is named:", required=True, type = str)
    parser.add_argument("-r2", "--read2", help="Your read 2 file of type fastq is named:", required=True, type = str)
    parser.add_argument("-q", "--qscore", help="Your quality score thershold:", required=True, type = int)
    parser.add_argument("-r3", "--read3", help="Your read 3 file of type fastq is named:", required=True, type = str)
    parser.add_argument("-r4", "--read4", help="Your read 4 file of type fastq is named:", required=True, type = str)
    parser.add_argument("-l", "--length", help="The length of each read is", required=True, type = int)
    parser.add_argument("-p", "--plot", help="The file you want to write your plot to is", required=True, type = str)
    return parser.parse_args()

args = get_args()

# print(f'barcode file name: {args.barcode_list_file}')

i_library_forward: set = set() #a set to hold the index library
i_library_rc: set = set()

def rc(fw_seq: str) -> str:
    '''A function to take the reverse complement of a sequence.'''
    rv_seq = ''
    rv_seq = [fw_seq[len(fw_seq)-1-k] for k, char in enumerate(fw_seq)] #reversing the sequence
    rc_seq: str = ''
    for base in rv_seq: #Taking complement of the strand
        if base == 'A':
            base = 'T'
        elif base == 'T':
            base = 'A'
        elif base == 'C':
            base = 'G'
        elif base == 'G':
            base = 'C'
        rc_seq += base
    return rc_seq


#### Main

with open(args.barcode_list_file, "r") as bf:
    for i, line in enumerate(bf):
        line = line.strip()
        if i > 0:
            seq = line.split()[4]
            i_library_forward.add(line.split()[4]) #Adds the indexes to a set
            i_library_rc.add(rc(seq))

i1: str = '' #index 1
i2: str = '' #index 2

c = 0 #file line counter

header1: str = '' #header line of first file
read1: str = '' #sequence line of first file (biological read 2)
qscore1: str = '' #qscore of first file

header2: str = '' #header line of second file
qscore2: str = '' #qscore of second file

header3: str = '' #header line of third file
qscore3: str = '' #qscore of third file

header4: str = '' #header line of fourth file
read2: str = '' #sequence line of fourth file (biological read 2)
qscore4: str = '' #qscore of fourth file

qscore_threshold: int = args.qscore #Quality score threshold. No index with a mean qscore less than this will be saved.


check = 0

R1_matched_files_dict: dict = {}
R2_matched_files_dict: dict = {}
#make a dictionary whos keys are the 24 indexes and the value is the whole file handle.

for indexes in i_library_forward: #Creates and opens 48 files for writing.
    #To call the files for writing later, use: "R1_matched_files_dict[indexes].write()"
    R1_matched_files_dict[indexes] = open(f'Final_dir20_real/R1_{indexes}.fastq', "w")
    R2_matched_files_dict[indexes] = open(f'Final_dir20_real/R2_{indexes}.fastq', "w")


eof: bool = False
#Need a count of each possible hopping combination.
#Use a dictionary where the keys are tuples of length 2 containing the two hopped indexes and the value is a count which increments.

fw_hopped_count = 0
rv_hopped_count = 0

fw_unknown_count = 0
rv_unknown_count = 0


matched_count_dict: dict = {}
matched_indexes_set: set = set()
matched_percent_list: list = []

hopped_dict: dict = {}
hopped_indexes_set: set = set()
hopped_percent_list: list = []

unknown_count: int = 0
hopped_count: int = 0
matched_count: int = 0

index_count_dict: dict = {}

total_read_counter: int = -1

index_percent_list: list = []
with gzip.open(args.read1, "rt") as r1f, gzip.open(args.read2, "rt") as r2f, gzip.open(args.read3, "rt") as r3f, gzip.open(args.read4, "rt") as r4f, open("Final_dir20_real/fw_hopped.fastq", "w") as fw_h, open("Final_dir20_real/rv_hopped.fastq", "w") as rv_h, open("Final_dir20_real/fw_unknown.fastq", "w") as fw_u, open("Final_dir20_real/rv_unknown.fastq", "w") as rv_u:
    while True:
        written: bool = False
        if eof == False:
            for c1 in range(4):
                l1 = r1f.readline().strip()
                if not l1:
                    eof = True
                    break
                if (c1 == 0):
                    # print("This is the header line")
                    header1 = l1
                elif (c1 == 1):
                    # print("This is the biological read 1 sequence line")
                    read1 = l1
                elif (c1 == 3):
                    # print("This is the quality score line")
                    qscore1 = l1
        if eof == False:
            for c2 in range(4):
                l2 = r2f.readline().strip()
                if not l2:
                    eof = True
                    break
                if (c2 == 0):
                    # print("This is the header line")
                    header2 = l2
                elif (c2 == 1):
                    # print("This is the index 1 line")
                    i1 = l2
                elif (c2 == 3):
                    # print("This is the quality score line")
                    qscore2 = l2
            
        if eof == False:
            for c3 in range(4):
                l3 = r3f.readline().strip()
                if not l3:
                    eof = True
                    break
                if (c3 == 0):
                    # print("This is the header line")
                    header3 = l3
                if (c3 == 1):
                    # print("This is the index 2 line")
                    i2 = rc(l3)
                if (c3 == 3):
                    # print("This is the quality score line")
                    qscore3 = l3
            
        if eof == False:
            for c4 in range(4):
                l4 = r4f.readline().strip()
                if not l4:
                    eof = True
                    break
                if (c4 == 0):
                    # print("This is the header line")
                    header4 = l4
                if (c4 == 1):
                    # print("This is biological read 2 line")
                    read2 = l4
                if (c4 == 3):
                    # print("This is the quality score line")
                    qscore4 = l4
            
        # print(f'bio read 1: {read1}\nindex 1: {i1}\n index 2: {i2}\n bio read 2: {read2}')
        if i1 not in i_library_forward or i2 not in i_library_forward or bioinfo.qual_score(i1) < qscore_threshold or bioinfo.qual_score(i2) < qscore_threshold and not eof: #separated the two if statements to filter out N's before comparing qual score to threshold
            #write to appropriate unknown files
           
            fw_u.write(f'{header1} {i1}-{i2}\n{read1}\n+\n{qscore1}\n')
            rv_u.write(f'{header4} {i1}-{i2}\n{read2}\n+\n{qscore4}\n')
            unknown_count += 1
        else:
            if i1 == i2: #Checking to see if i1 and i2 match.
                #Write to appropriate matched files
                R1_matched_files_dict[i1].write(f'{header1} {i1}-{i2}\n{read1}\n+\n{qscore1}\n')
                R2_matched_files_dict[i2].write(f'{header4} {i1}-{i2}\n{read2}\n+\n{qscore4}\n')
               
                matched_index_tuple: tuple = (i1, i2)
                if matched_index_tuple not in matched_count_dict:
                    matched_count_dict[matched_index_tuple] = 1
                else:
                    matched_count_dict[matched_index_tuple] += 1
                matched_indexes_set.add(matched_index_tuple)
                matched_count += 1

            else: #if they don't match then they must have hopped.
                #Write to appropriate hopped files
          
                fw_h.write(f'{header1} {i1}-{i2}\n{read1}\n+\n{qscore1}\n')
                rv_h.write(f'{header4} {i1}-{i2}\n{read2}\n+\n{qscore4}\n')
                
            
                hopped_count += 1

                hopped_index_tuple: tuple = (i1, i2)
                if hopped_index_tuple not in hopped_dict:
                    hopped_dict[hopped_index_tuple] = 1
                else:
                    hopped_dict[hopped_index_tuple] += 1
                hopped_indexes_set.add(hopped_index_tuple)

        total_read_counter += 1

        if i1 not in index_count_dict:
            index_count_dict[i1] = 1
        else:
            index_count_dict[i1] += 1
        if i2 not in index_count_dict:
            index_count_dict[i2] = 1
        else:
            index_count_dict[i2] += 1

        if eof == True:
            break

#Calculating and writing the statistics to the stats file
with open("Final_dir20_real/Demultiplexing_stats_file", "w") as dsf:
    dsf.write(f'Index-pair\tCount\tPercent of matched reads by index-pair\n')
    for i, hopped_i in enumerate(hopped_indexes_set):
        hopped_percent_list.append((hopped_dict[hopped_i]/total_read_counter)*100)
        for z, indv_hopped_index in enumerate(hopped_i):
            dsf.write(indv_hopped_index)
            if z < 1:
                dsf.write("-")
        dsf.write(f': {hopped_dict[hopped_i]}\t{hopped_percent_list[i]}\n')
    for i, matched_i in enumerate(matched_indexes_set):
        matched_percent_list.append((matched_count_dict[matched_i]/total_read_counter)*100)
        
        for z, matched_index in enumerate(matched_i):
            dsf.write(matched_index)
            if z < 1:
                dsf.write("-")
        dsf.write(f'{matched_count_dict[matched_i]}\t{matched_percent_list[i]}')
    dsf.write(f'Percent of total reads that matched: {(matched_count/total_read_counter)*100}\n')
    dsf.write(f'Percent of total reads that hopped: {(hopped_count/total_read_counter)*100}\n')
    dsf.write(f'Percent of total unknown reads: {(unknown_count/total_read_counter)*100}\n')
    dsf.write(f'Total number of reads: {(total_read_counter)}\n')

#Making a bar plot of frequency of each index
index_count_list = list(index_count_dict.values())

for index_count_value in index_count_list:
    index_percent_list.append((index_count_value/total_read_counter)*100) # y-values in plot

index_list = list(index_count_dict.keys()) #x-values in plot

fig, ax = plt.subplots(figsize=(14, 6))

ax.bar(index_list, index_percent_list, color = "purple")
ax.set_title('Frequency of each index across all reads', fontsize = 16)
ax.set_xlabel('Index', fontsize = 14)
ax.set_ylabel('Percent of reads with that index', fontsize = 14)

plt.tight_layout()
# plt.show()

plt.savefig(args.plot)
plt.close()



for indexes in i_library_forward: #Creates and opens 48 files for writing.
    #To call the files for writing later, use: "R1_matched_files_dict[indexes].write()"
    R1_matched_files_dict[indexes].close()
    R2_matched_files_dict[indexes].close()
#NOTE you cannot explicitly check for eof using a for loop.
#NOTE 