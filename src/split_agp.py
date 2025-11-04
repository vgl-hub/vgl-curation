#!/usr/bin/env python3

import textwrap
import os
import argparse
import csv
import re
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
from binascii import hexlify
import pandas as pd
import sys
import gzip
import csv



def main():
    parser = argparse.ArgumentParser(
                        prog='split_agp',
                        description='Correct AGP for sequence lengths, split the agp per haplotype, assign unlocs and remove duplicated haplotigs  ',
                        usage='split_agp -f assembly.fasta -a curated_agp.agp',
                        formatter_class=argparse.RawTextHelpFormatter,
                        epilog=textwrap.dedent('''
                                            Outputs:
                                            - corrected.agp: Curated AGP file corrected for sequence lengths.
                                            - Hap_1/hap1.agp: AGP file for Haplotype 1.
                                            - Hap_2/hap2.agp: AGP file for Haplotype 2.
                                            - Hap_1/hap.unlocs.no_hapdups.agp: AGP file for Haplotype 1 with assigned unlocs and purged of duplicated haplotigs.
                                            - Hap_2/hap.unlocs.no_hapdups.agp: AGP file for Haplotype 2 with assigned unlocs and purged of duplicated haplotigs.
                                            - Hap_1/haplotigs.agp: AGP file containing the haplotig duplications removed from Haplotype 1.
                                            - Hap_2/haplotigs.agp: AGP file containing the haplotig duplications removed from Haplotype 2.
                                            '''))

    parser.add_argument('-a', '--agp', dest="agp", required=True, help='Path to the curated AGP file')
    parser.add_argument('-f', '--fasta', dest="fasta",required=True, help='Path to the assembly fasta file')
    parser.add_argument('-o', '--output', dest="output_dir", default="./", required=False, help='Path to the output directory (Default: current directory)') 
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)


    AGPcorrect(args.fasta, args.agp, args.output_dir+'corrected.agp')

    hap_split(args.output_dir+'corrected.agp',args.output_dir+'Hap_1/hap1.agp',args.output_dir+'Hap_2/hap2.agp')

    print("\nRun unloc on Haplotype 1\n")
    unloc(args.output_dir+'Hap_1/hap1.agp',args.output_dir+'Hap_1/')
    print("\nRun unloc on Haplotype 2\n")
    unloc(args.output_dir+'Hap_2/hap2.agp',args.output_dir+'Hap_2/')


if __name__ == "__main__":
    main()
        



def Open(file_name):
    with open(file_name, "rb") as f:
        isgzip = hexlify(f.read(2)) == b"1f8b"
    return gzip.open(file_name, "rt") if isgzip else open(file_name, "r")


def AGPcorrect(fasta_file, agp_file, corrected_agp):

    #if len(sys.argv) == 1 or (sys.argv[1] in ("-h", "--help")):
    #    print(
    #        "Usage: AGPCorrect ref.fa(.gz) scaffs.agp > corrected_scaffs.agp",
    #        file=sys.stderr,
    #    )
    #    sys.exit(0)

    #if len(sys.argv) != 3:
    #    sys.exit("Usage: AGPCorrect ref.fa(.gz) scaffs.agp > corrected_scaffs.agp")

    with Open(fasta_file) as f:
        seqs = {seq.id: len(seq) for seq in SeqIO.parse(f, "fasta")} ## Generating a list of all sequences and their lengths from the fasta
        ## All expected sequences are present. 

    # print(
    #     f"Read fasta, {len(seqs)} sequences",
    #     *(f"{s}: {n} bp" for s, n in seqs.items()),
    #     "\n",
    #     file=sys.stderr,
    #     sep="\n",
    # )

    seen = {}
    with open(agp_file, "r") as f: ## Building a dictionary for scaffolds and their lengths based on AGP output
        for line in f:
            if not line.startswith("#"):
                line = line.split("\t")
                if line[4] == "W":
                    seen[line[5]] = max(seen.setdefault(line[5],0), int(line[7]))
    ## everything printed here is expected


    #stdout_file=sys.stdout
    write_out= open(corrected_agp, 'w')
                    
    with open(agp_file, "r") as f: ## Opening the AGP again
        curr_scaff = None 
        maxn = 1
        for line in f:
            line = line[:-1]
            if not line.startswith("#"):
                line = line.split("\t")
                if curr_scaff != line[0]: ## checking whether the current scaffold name matches the scaffold name in the previous line
                ## If the current scaffold name does not match the previous, then it print ths scaffold name and the "correction" (which I assume is the BP difference between the )
                    # if curr_scaff:
                    #     print(f"{curr_scaff}: {correct} bp correction", file=sys.stderr)
                    curr_scaff = line[0]
                    # print (line, curr_scaff, "<- AGP")
                    correct = 0 

                line[1] = str(int(line[1]) + correct) ## Printing corrected start position I believe

                if line[4] == "W" and ((this_l := int(line[7])) == seen[line[5]]):
                    correct += (acc_l := seqs[line[5]]) - this_l ## everything seems to be fine here, the sizes haven't been corrected 
                
                
                    ##seems like this is where we are getting negative values from - but it could be that the right "correction isn't matched with the right scaffold "
                    # print (seqs[line[5]], line[5]) 
                    # print (this_l, correct, '\n')
                    ## Seqs is a dictionary of scaffolds and their lengths generated from the original fasta file 
                    ## So line[5] is the name of the current scaffold and seqs[line[5]] goes to the dict and produces the true length of the scaffold 

                    if int(line[6]) >= acc_l:  ##line[6] is the start position of each segment - seems to be checking whether the start position is greater than the actual length of the scaffold?
                        sys.exit(
                            "Error with line: {}\n{} > {}".format(
                                "\t".join(line), line[6], acc_l
                            )
                        )

                    line[7] = str(acc_l)

                line[2] = str(int(line[2]) + correct)

                print("\t".join(line), file=write_out)
                maxn = max(maxn, int(line[0].split("_")[-1]))
            else:
                if line.startswith("# DESCRIPTION"):
                    line += "\tModified by PretextView_AGPCorrect"

        # if curr_scaff:
        #     print(f"{curr_scaff}: {correct} bp correction", file=sys.stderr)


    maxn += 1
    for k, (s, n) in enumerate((s, n) for s, n in seqs.items()):
        if s not in set(seen.keys()):
            print(f"Scaffold_{maxn+k}\t1\t{n}\t1\tW\t{s}\t1\t{n}\t+", file=write_out)

    write_out.close()



def hap_split(corrected_agp,path_agp1,path_agp2):


    corr=corrected_agp

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
            if 'H1' in line[5] or 'hap1' in line[5] :
                H1_lines.append(line)
            elif 'H2' in line[5] or 'hap2' in line[5]:
                H2_lines.append(line)
            else:
                print("Warning: Scaffold "+line[5]+" not assigned to any haplotype. It will be skipped. Verify your scaffolds have haplotype markers type 'H1' or 'hap1'.", file=sys.stderr)

    #Create output directories
    if '/' in path_agp1:
        tmp_list1=path_agp1.split('/')
        output_dir1="/".join(tmp_list1[:-1])
    os.makedirs(output_dir1, exist_ok=True)

    if '/' in path_agp2:
        tmp_list2=path_agp2.split('/')
        output_dir2="/".join(tmp_list2[:-1])
    os.makedirs(output_dir2, exist_ok=True)


    with open (path_agp1,'w',newline='\n') as file1:
        writer=csv.writer(file1,delimiter='\t')
        writer.writerows(header)
        writer.writerows(H1_lines)
    file1.close()

    with open (path_agp2, 'w', newline='\n') as file2:
        writer=csv.writer(file2,delimiter='\t')
        writer.writerows(header)
        writer.writerows(H2_lines)
    file2.close()




def unloc(hap_agp,output_dir):

    outdir=output_dir
    hap=hap_agp

    os.makedirs(outdir, exist_ok=True)

    header=[]
    agp_lines=[]
    with open(hap) as file:
        agp = csv.reader(file,delimiter='\t')
        for line in agp:
            if "#" in line[0]:
                header.append(line)
            else: 
                agp_lines.append(line)
    file.close()

    maxlen=max([len(entry) for entry in agp_lines])

    if maxlen < 11:
        ## Checking for presence of metadata tags - if no unlocs or haplotigs present, the script just generates a replicate agp and exits.
        print ("No metadata tags used. Are you sure there are no unlocs, haplotigs or sex chromosomes to label?")
        with open (outdir+'/hap.unlocs.no_hapdups.agp','w',newline='\n') as f:
            writer=csv.writer(f,delimiter='\t')
            writer.writerows(header)
            writer.writerows(agp_lines)
        f.close() 
        exit()
        
        
    agp_df=pd.DataFrame(agp_lines)

    agp_df=agp_df.rename(columns={0:'chr',1:'chr_start', 2:'chr_end', 3:'#_scaffs', 4:'W', 5:'scaff', 6:'scaff_start', 7:'scaff_end', 8:'ori', 9:'painted', 10:'tag'})



    unlocs=(agp_df.index[agp_df['tag']=='Unloc']).to_list()

    scaffs_with_unlocs=[]
    unloc_num=1
    for index in unlocs: ##This assumes unlocs are placed at the end of scaffolds 
        scaff=(agp_df.iloc[index].to_list())[0]
        agp_df.loc[index,'chr_start']=1
        agp_df.loc[index,'chr_end']=agp_df.loc[index,'scaff_end']
        if scaff in scaffs_with_unlocs:
            unloc_num+=1
            agp_df.loc[index,'chr']=agp_df.loc[index,'chr']+"_unloc_"+str(unloc_num)
        else:
            unloc_num=1
            agp_df.loc[index,'chr']=agp_df.loc[index,'chr']+"_unloc_"+str(unloc_num)
            scaffs_with_unlocs.append(scaff)

    haplotigs=(agp_df.index[agp_df['tag']=='Haplotig'])
    agp_df_mod=agp_df.drop(haplotigs)
    agp_list=agp_df_mod.values.tolist()

    num_lines=len(agp_list)-1
    line_num=0
    prox_lig_lines=[]
    while line_num < num_lines:
        current_line=agp_list[line_num]
        current_scaff=current_line[0]
        prev_line=agp_list[line_num-1]
        next_line=agp_list[line_num+1]
        if current_line[10]=='Unloc' and prev_line[8]=='proximity_ligation':
            prox_lig_lines.append(prev_line)
        elif current_line[10]=='Unloc' and next_line[8]=='proximity_ligation':
            prox_lig_lines.append(next_line)
        elif current_line[8]=='proximity_ligation' and prev_line[8]=='proximity_ligation':
            prox_lig_lines.append(current_line)
            prox_lig_lines.append(prev_line)
        elif current_line[0]!=next_line[0] and current_line[8]=='proximity_ligation':
            prox_lig_lines.append(current_line)
        line_num+=1


    final_list=[]
    for line in agp_list:
        if line in prox_lig_lines:
            print("Gap line removed: "+ "\t".join(str(n) for n in line)+"\n")
        else:
            final_list.append(line)

    haplotigs_list=[agp_df.iloc[ind].tolist() for ind in haplotigs]


    with open ((outdir+'/hap.unlocs.no_hapdups.agp'),'w',newline='\n') as f:
        writer=csv.writer(f,delimiter='\t')
        writer.writerows(header)
        writer.writerows(final_list)
    f.close()

    with open ((outdir+'/haplotigs.agp'),'w',newline='\n') as h:
        writer=csv.writer(h,delimiter='\t')
        writer.writerows(haplotigs_list)
    h.close()




