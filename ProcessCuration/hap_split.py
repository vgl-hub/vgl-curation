#!/usr/bin/env python3

## for splitting the haplotypes
## After AGPcorrect but prior to unloc assignment 

import csv 
import pandas as pd
import argparse
import textwrap
import os

parser = argparse.ArgumentParser(
                    prog='hap_split.py',
                    description='Separate the AGP file for two haplotypes into one AGP file for each haplotype. ',
                    usage='hap_split.py -1 Hap_1/hap.agp -2 Hap_2/hap.agp -a curated_agp.agp ',
                    formatter_class=argparse.RawTextHelpFormatter,
                    epilog=textwrap.dedent('''
                                           Outputs: 
                                           - {hap1_directory}/hap.agp: AGP file for Haplotype 1
                                           - {hap2_directory}/hap.agp: AGP file for Haplotype 2
                                           '''))
 
parser.add_argument('-a', '--agp', dest="agp", required=True, help='Path to the curated AGP file')
parser.add_argument('-1', '--hap1', dest="hap1",required=True, help='Path to the AGP file for haplotype 1')  
parser.add_argument('-2', '--hap2', dest="hap2", required=True, help='Path to the AGP file for haplotype 2') 
args = parser.parse_args()


corr=args.agp

header=[]
agp_lines=[]
with open(corr) as file:
    agp = csv.reader(file,delimiter='\t')
    for line in agp:
        if "#" in line[0]:
            header.append(line)
        else: 
            agp_lines.append(line)
file.close()


current_hap=''
current_scaff=''
H1_lines=[]
H2_lines=[]
for line in agp_lines:
    
    if 'proximity_ligation' in line or 'Painted' in line:
        if 'Hap_1' in line:
            H1_lines.append(line)
            current_hap='Hap_1'
        elif 'Hap_2' in line: 
            H2_lines.append(line)
            current_hap='Hap_2'
        elif current_hap=='Hap_1':
            H1_lines.append(line)
        elif current_hap=='Hap_2':
            H2_lines.append(line)
    else:
        if 'H1' in line[5]:
            H1_lines.append(line)
        elif 'H2' in line[5]:
            H2_lines.append(line)

#Create output directories
if '/' in args.hap1:
    tmp_list1=args.hap1.split('/')
    output_dir1="/".join(tmp_list1[:-1])
os.makedirs(output_dir1, exist_ok=True)

if '/' in args.hap2:
    tmp_list2=args.hap2.split('/')
    output_dir2="/".join(tmp_list2[:-1])
os.makedirs(output_dir2, exist_ok=True)


with open (args.hap1,'w',newline='\n') as file1:
    writer=csv.writer(file1,delimiter='\t')
    writer.writerows(header)
    writer.writerows(H1_lines)
file1.close()

with open (args.hap2, 'w', newline='\n') as file2:
    writer=csv.writer(file2,delimiter='\t')
    writer.writerows(header)
    writer.writerows(H2_lines)
file2.close()
