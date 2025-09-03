#!/usr/bin/env python3

import pandas as pd 
import argparse
from io import StringIO
from natsort import natsorted
import textwrap
import os
import re
import sys


def main():

    parser = argparse.ArgumentParser(
                        prog='filter_mashmap_with_tagged_pairs',
                        description='Filter Mashmap output to keep only Scaffolds paired with the tags Hap_1 and Hap_2',
                        usage='filter_mashmap_with_tagged_pairs -1 Hap_1/inter_chr.tsv -2 Hap_2/inter_chr.tsv -q Hap_2 -r Hap_1 -agp curated_agp_with_micro_tags.agp -m mashmap.out -o results/ ',
                        formatter_class=argparse.RawTextHelpFormatter,
                        epilog=textwrap.dedent('''
                                            Outputs: 
                                            - {out_dir}/rvcp.sak: SAK file for gfastats reversing of hap2 sequences
                                            - {out_dir}/orientation.tsv:  Table with the orientation of Hap1 vs Hap2 scaffolds (The sign is the most represented in the mashmap alignments between two sequences)
                                            - {out_dir}/hap2.vs.hap1.tsv: The Mapping between hap2 and hap1 names. 
                                            '''))
    
    parser.add_argument('-1', '--hap1', dest="hap1",required=True, help='Path to the chromosome assignment file for  Haplotype 1 (inter_chr.tsv) ')  
    parser.add_argument('-2', '--hap2', dest="hap2", required=True, help='Path to the chromosome assignment file for Haplotype 2 (inter_chr.tsv) ') 
    parser.add_argument('-q', '--query', dest="query", default="Hap_2", help='Haplotype use as Query for MashMap: Hap_1 or Hap_2 (Default Hap_2)')  
    parser.add_argument('-r', '--reference', dest="reference", default="Hap_1", help='Haplotype use as reference for MashMap: Hap_1 or Hap_2 (Default Hap_1)')  
    parser.add_argument('-a', '--agp', dest="agp", required=True, help='Path to the curated AGP file (with the tags Hap_1 and Hap_2)')       
    parser.add_argument('-m', '--mashmap', dest="mashmap", required=True, help='Path to the curated Mashmap output')  
    parser.add_argument('-o', '--out_dir', dest="out_dir", required=True, help='Path to output Prefix')  
    args = parser.parse_args()


    if '/' in args.out_dir:
        tmp_list1=args.out_dir.split('/')
        output_dir1="/".join(tmp_list1[:-1])
    os.makedirs(output_dir1, exist_ok=True)


    if args.query!="Hap_1" and args.query!="Hap_2":
        raise SystemExit("The query parameter should be Hap_1 or Hap_2")
    if args.reference!="Hap_1" and args.reference!="Hap_2":
        raise SystemExit("The reference parameter should be Hap_1 or Hap_2")

    try:
        dico_supers_hap1= pd.read_csv(args.hap1, sep=r'\s+', header=None)
        dico_supers_hap2= pd.read_csv(args.hap2, sep=r'\s+', header=None)
    except FileNotFoundError:
        print(f"Error: The file inter_chr.tsv was not found.")

    try:
        mashmap_out=pd.read_csv(args.mashmap, sep=r'\s+', header=None)
    except FileNotFoundError:
        print(f"Error: The file '{args.mashmap}' was not found.")

    mashmap_out=mashmap_out.rename(columns={0:args.query, 1:'Query length',2:'Query 0-based start',3:'Query 0-based end',4:'Orientation',5:args.reference, 6:'Target length', 7:'Target 0-based start', 8:'Target 0-based end',12:'Mapping nucleotide identity'})

    filtered_mashmap=pd.DataFrame(columns=mashmap_out.columns)

    file_agp = args.agp
    filter_word = "Hap_"
    output_agp=""
    try:
        with open(file_agp, 'r') as file:
            for line in file:
                if filter_word in line:
                    tmp_line_split=line.split(r'\s+')
                    line=re.sub(r"Painted\s+([M])",r"Painted\tNA\tM",line)
                    line=re.sub(r"Painted\s+([^HN])",r"Painted_\1",line)
                    output_agp=output_agp+line
    except FileNotFoundError:
        print(f"Error: The file '{file_name}' was not found.")

    filtered_agp = pd.read_table(StringIO(output_agp),sep=r'\s+',header=None)
    if len(filtered_agp[filtered_agp[9]=="Painted_Y"])!=0 or len(filtered_agp[filtered_agp[9]=="Painted_W"])!=0:
        filtered_agp=filtered_agp[~filtered_agp[9].str.contains('Painted_')]

    filtered_agp=filtered_agp[[0,5,8,9,10]]
    duplications=filtered_agp[filtered_agp.duplicated()]
    if len(duplications)>0:
        print("Error: The following scaffolds are duplicated in the agp: ", file=sys.stderr)
        print(duplications.to_string(index=False, header=False), file=sys.stderr) 
        raise SystemExit("Please verify your curation.")

    filtered_agp=filtered_agp.reset_index(drop=True)

    dico_micro=pd.DataFrame(columns=['Hap_1', 'Hap_2'])

    even_value = None
    even_hap=None

    for index, row in filtered_agp.iterrows():  
        if index % 2 == 0:
            even_value=row[0] 
            even_hap=row[10]
        else:
            index_new=len(dico_micro)
            tmp=pd.DataFrame(columns=['Hap_1', 'Hap_2'])
            tmp.loc[0,row[10]]=row[0]
            tmp.loc[0,even_hap]=even_value
            dico_micro.loc[len(dico_micro)] = tmp.loc[0]

    Paired = pd.DataFrame(columns=['Hap_1', 'Hap_2'])

    try:
        for index, row in dico_micro.iterrows():
            Scaffold1=row.iloc[0]
            Scaffold2=row.iloc[1]
            res_search1=dico_supers_hap1[dico_supers_hap1[0]==Scaffold1]
            res_search2=dico_supers_hap2[dico_supers_hap2[0]==Scaffold2]
            if len(res_search1)==0 or len(res_search2)==0:                  
                raise SystemExit("Error: Names from the curated agp were not founds in the hap.unlocs.no_hapdups.agp file. Verify that the names are consistants and all painted scaffolds have the tags Hap_1 or Hap_2")
            super1=dico_supers_hap1[dico_supers_hap1[0]==Scaffold1].iloc[0,1]
            super2=dico_supers_hap2[dico_supers_hap2[0]==Scaffold2].iloc[0,1]
            Paired.loc[index] = [super1, super2]
    except Exception as e:
        print(f"An error occurred: {e}")

    for index, row in mashmap_out.iterrows():
        if  ((Paired['Hap_1'] == row[[args.query,args.reference]].iloc[1]) & (Paired["Hap_2"] == row[[args.query,args.reference]].iloc[0])).any():
            filtered_mashmap.loc[len(filtered_mashmap)] = row

    orientations = filtered_mashmap[["Hap_1","Hap_2","Orientation"]]



    # Group by Reference and Query, then count Sign occurrences
    counts = orientations.groupby(["Hap_1","Hap_2"])["Orientation"].value_counts()

    # Get the sign with the highest count per group
    most_frequent_sign = counts.groupby(level=[0, 1]).idxmax().apply(lambda x: x[2])  # x is a tuple (Reference, query, sign)

    # Combine into a final dataframe
    result = most_frequent_sign.to_frame(name='Main Orientation')
    result=result.reset_index()

    # Merge with indicator to see where rows originate
    merged_df = Paired[['Hap_1','Hap_2']].merge(result[['Hap_1','Hap_2']], how='left', indicator=True)
    # Filter for rows not present in Mashmap
    rows_not_in_mashmap = merged_df[merged_df['_merge'] == 'left_only']
    if len(rows_not_in_mashmap)>0:
        rows_not_in_mashmap[['Hap_1','Hap_2']].to_csv(args.out_dir+"/missing.tsv", index=False, header=True, sep="\t")
        print("Warning: One or more pairs of scaffolds identified by curation are not found in mashmap. See "+args.out_dir+"missing.tsv" )
        rows_not_in_mashmap= rows_not_in_mashmap.copy()
        rows_not_in_mashmap.loc[:, "Main Orientation"] = "+"
        result = pd.concat([result, rows_not_in_mashmap[["Hap_1","Hap_2","Main Orientation"]]])


    result['Hap_1'] = pd.Categorical(result['Hap_1'], categories=natsorted(result['Hap_1'].unique()), ordered=True)
    result = result.sort_values(by='Hap_1', ignore_index=True)

    results_with_unloc=result.copy()

    for i, row in result.iterrows():
        super_name2=row['Hap_2']
        super_name1=row['Hap_1']
        unloc_search=dico_supers_hap2[dico_supers_hap2[1].str.contains(super_name2+"_unloc", case=False, na=False)]
        if len(unloc_search)!=0 : 
            for j, unloc in unloc_search.iterrows():
                unloc_old_name=unloc[1]
                unloc_new_name=unloc_old_name.replace(super_name2, super_name1)
                results_with_unloc.loc[len(results_with_unloc)] = [unloc_new_name, unloc_old_name, "+"]

        
    results_with_unloc.to_csv(args.out_dir+"/orientation.tsv", index=True, header=True, sep="\t")
    
    renaming=results_with_unloc[results_with_unloc['Hap_1']!=results_with_unloc['Hap_2']][['Hap_1','Hap_2']]
    renaming['Hap_2'] = renaming['Hap_2'].astype(str) + "_oldname"
    renaming['Command']="RENAME"


    results_with_unloc['Hap_2'] = results_with_unloc['Hap_2'].astype(str) + "_oldname"

    to_reverse = pd.DataFrame(columns=['Command', 'Hap_2'])
    to_reverse['Hap_2']=results_with_unloc[results_with_unloc['Main Orientation']=='-']["Hap_2"]
    to_reverse['Command']="RVCP"


    #to_reverse.to_csv(args.out_dir+"/rvcp.sak", index=False, header=False, sep="\t")

    #renaming=renaming.rename(columns={'Hap_1':'New_name','Hap_2':'Old_name'})

    #renaming[['Command','Hap_2','Hap_1']].to_csv(args.out_dir+"/hap2.vs.hap1.tsv", index=False, header=false, sep="\t")

    to_reverse['Hap_1']=''

    actions= pd.concat([to_reverse, renaming])
    actions[['Command','Hap_2','Hap_1']].to_csv(args.out_dir+"/reversing_renaming.sak", index=False, header=False, sep="\t")    



if __name__ == "__main__":
    main()
        