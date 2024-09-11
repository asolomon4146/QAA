#!/usr/bin/env python

import bioinfo
import argparse
import matplotlib.pyplot as plt
import gzip
import numpy as np

def get_args():
    parser = argparse.ArgumentParser(description="A program to extract statistics from demultiplexed reads")
    parser.add_argument("-i", "--input", help="Your input data file name is", required=True, type = str)
    parser.add_argument("-l", "--length", help="The length of each read is", required=True, type = int)
    parser.add_argument("-p", "--plot", help="The file you want to write your plot to is", required=True, type = str)
    return parser.parse_args()

args = get_args()
read_length: int = int(args.length)

input_file = args.input

def init_list(lst: list, value: float=0.0) -> list:
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with the number found in args.length of 0.0.'''
    
    for i in range(0,read_length):
        lst.append(value)
    return lst

my_list: list = []
my_list = init_list(my_list)

qscore_array = np.zeros((101, 1000000))

def populate_list(input_file: str, qscore_array) -> tuple[list, int]:
    """Creates an empty list, opens the fastq file, reads it line by line,
    parses out the quality scores, converts them to numbers, sums the Qscores in each base, returns the list with a counter"""
    i = 0
    seq: str = ""
    seq_count: int = 0
    sum_Q = init_list([], 0.0)
    with gzip.open(input_file, "rt") as f:
        for line in f:
            line = line.strip()
            if (i+1) %4 == 2:
                seq = line
            if "N" not in seq: #Does not evaluate any indexes which contain an N.
                if (i+1) %4 == 0:
                    for c, base in enumerate(line):
                        sum_Q[c]+=(bioinfo.convert_phred(base))
                        if seq_count < 41000000 and seq_count >= 40000000:
                            qscore_array[c, seq_count-40000000]=(bioinfo.convert_phred(base))
                    seq_count += 1

            i+=1
    return sum_Q, seq_count

my_list, num_lines = populate_list(input_file, qscore_array)
print("\n\nSeq count:", num_lines, "\n\n")

my_list = [j/(num_lines) for j in my_list]
num_reads = num_lines

def calc_median(qscores):
    func_med = np.zeros(101)
    for i in range(len(qscores)):
        func_med[i] = np.median(qscores[i])
    return(func_med)

median = calc_median(qscore_array)

def calc_stdev(qscores):
    func_stdev = np.zeros(101)
    for i in range(len(qscores)):
        func_stdev[i] = np.std(qscores[i], ddof = 1)
    return(func_stdev)

stdev = calc_stdev(qscore_array)

print("# Base Pair" + "\t" + "Mean Quality Score" + "\t" "Median Quality Score" + "\t" "Standard Deviation")
for k, mqs in enumerate(my_list):
    print(f'{k}\t{mqs}\t{median}\t{stdev}')


# def calc_stdev(qscores, length: int):
#     func_stdev = np.zeros(length)
#     for i in range(len(qscores)):
#         func_stdev[i] = np.std(qscores[i], ddof = 1)
#     return(func_stdev)

# stdev = calc_stdev(my_list, read_length)


#Plotting with matplotlib:
print("Object type of read_length variable:", type(read_length), "\nAnd the value of read_length is:", read_length)
fig, ax = plt.subplots(figsize=(14, 6))

positions = list(range(1, (read_length+1)))

ax.plot(positions, my_list, ".-", color = "blue", linewidth = 1)
# ax.errorbar(positions, my_list, yerr = stdev, color='blue', linewidth=1)
ax.set_title('Average Quality Score by Nucleotide Position', fontsize = 16)
ax.set_xlabel('Nucleotide Position', fontsize = 14)
ax.set_ylabel('Mean Quality Score', fontsize = 14)

plt.tight_layout()
# plt.show()

plt.savefig(args.plot)
# plt.close()

# fig, ax = plt.subplots(figsize=(14, 6))

# ax.errorbar(positions, median, yerr = stdev, color='red', linewidth=1)
# ax.set_title('Median Quality Score by Nucleotide Position', fontsize = 16)
# ax.set_xlabel('Nucleotide Position', fontsize = 14)
# ax.set_ylabel('Mean Quality Score', fontsize = 14)
# plt.tight_layout()
# plt.show()

