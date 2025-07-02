## for splitting the haplotypes
## After AGPcorrect but prior to unloc assignment 

import csv 
import pandas as pd
import argparse
import textwrap

parser = argparse.ArgumentParser(
                    prog='hap_split.py',
                    description='Separate the AGP file for two haplotypes into one AGP file for each haplotype. ',
                    usage='hap_split.py -1 Hap_1/ -2 Hap_2/ -a curated_agp.agp ',
                    formatter_class=argparse.RawTextHelpFormatter,
                    epilog=textwrap.dedent('''
                                           Outputs: 
                                           - {hap1_prefix}hap.agp: AGP file for Haplotype 1
                                           - {hap2_prefix}hap.agp: AGP file for Haplotype 2
                                           '''))
 
parser.add_argument('-a', '--agp', dest="agp", required=True, help='Path to the curated AGP file')
parser.add_argument('-1', '--hap1', dest="hap1",required=True, help='Output prefix for haplotype 1')  
parser.add_argument('-2', '--hap2', dest="hap2", required=True, help='Output prefix for haplotype 2') 
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


with open (args.hap1+'hap.agp','w',newline='\n') as file1:
    writer=csv.writer(file1,delimiter='\t')
    writer.writerows(header)
    writer.writerows(H1_lines)
file1.close()

with open (args.hap2+'hap.agp', 'w', newline='\n') as file2:
    writer=csv.writer(file2,delimiter='\t')
    writer.writerows(header)
    writer.writerows(H2_lines)
file2.close()
