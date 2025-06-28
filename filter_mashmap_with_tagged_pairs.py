import pandas as pd 
import argparse
from io import StringIO


parser = argparse.ArgumentParser(
                    prog='filter_mashmap_with_tagged_pairs',
                    description='Filter Mashmap output to keep only Scaffolds paired with the tags Micro1 and Micro2',
                    usage='filter_mashmap_with_tagged_pairs.py -1 Hap_1/ -2 Hap_2/ -m mashmap.out -o filtered_mashmap.out')
 
parser.add_argument('-1', '--hap1', dest="hap1",required=True, help='Folder containing the files for Haplotype 1 (hap.unlocs.no_hapdups.agp and inter_chr.tsv) ')  
parser.add_argument('-2', '--hap2', dest="hap2", required=True, help='Folder containing the files for Haplotype 2 (hap.unlocs.no_hapdups.agp and inter_chr.tsv) ') 
parser.add_argument('-a', '--agp', dest="agp", required=True, help='Path to the curated AGP file (containing the tags Micro1 and Micro2)')       
parser.add_argument('-m', '--mashmap', dest="mashmap", required=True, help='Path to the curated Mashmap output')  
parser.add_argument('-o', '--out', dest="out_file", required=True, help='Path output')  

args = parser.parse_args()

try:
    dico_scaffolds_hap1= pd.read_csv(args.hap1+"/hap.unlocs.no_hapdups.agp", sep="\t", header=None)
    dico_scaffolds_hap2= pd.read_csv(args.hap2+"/hap.unlocs.no_hapdups.agp", sep="\t", header=None)
except FileNotFoundError:
    print(f"Error: The file hap.unlocs.no_hapdups.agp was not found.")
    
try:
    dico_supers_hap1= pd.read_csv(args.hap1+"/inter_chr.tsv", sep="\t", header=None)
    dico_supers_hap2= pd.read_csv(args.hap2+"/inter_chr.tsv", sep="\t", header=None)
except FileNotFoundError:
    print(f"Error: The file inter_chr.tsv was not found.")

try:
    mashmap_out=pd.read_csv(args.mashmap, sep="\t", header=None)
except FileNotFoundError:
    print(f"Error: The file '{args.mashmap}' was not found.")

dico_scaffolds_hap2=dico_scaffolds_hap2[dico_scaffolds_hap2[11]=="Micro2"][[0,5]]
dico_scaffolds_hap1=dico_scaffolds_hap1[dico_scaffolds_hap1[11]=="Micro1"][[0,5]]

file_agp = args.agp
filter_word = "Micro"
output_agp=""
try:
    with open(file_agp, 'r') as file:
        for line in file:
            if filter_word in line:
                output_agp=output_agp+line
except FileNotFoundError:
    print(f"Error: The file '{file_name}' was not found.")
except Exception as e:
    print(f"An error occurred: {e}")

filtered_agp = pd.read_csv(StringIO(output_agp),sep='\t',header=None)

dico_micro=pd.DataFrame(columns=['Hap_1', 'Hap_2'])

even_value = None
even_hap=None

for index, row in filtered_agp.iterrows():  
    if index % 2 == 0:
        even_value=row.iloc[5]  # 0-based index; $6 is fields[5]
        even_hap=row.iloc[10]
    else:
        index_new=len(dico_micro)
        tmp=pd.DataFrame(columns=['Hap_1', 'Hap_2'])
        tmp.loc[0,row.iloc[10]]=row.iloc[5]
        tmp.loc[0,even_hap]=even_value
        dico_micro.loc[len(dico_micro)] = tmp.loc[0]


Paired = pd.DataFrame(columns=['Hap_1', 'Hap_2'])

try:
    for index, row in dico_micro.iterrows():
        Scaffold1=row.iloc[0]
        Scaffold2=row.iloc[1]
        res_search1=dico_scaffolds_hap1[dico_scaffolds_hap1[5]==Scaffold1]
        res_search2=dico_scaffolds_hap2[dico_scaffolds_hap2[5]==Scaffold2]
        if len(res_search1)==0 or len(res_search2)==0:                  
            raise ValueError("Names from the curated agp were not founds in the hap.unlocs.no_hapdups.agp file. Verify the names are consistants.")
        curated1=str(res_search1.iloc[0,0])
        curated2=str(res_search2.iloc[0,0])
        super1=dico_supers_hap1[dico_supers_hap1[0]==curated1].iloc[0,1]
        super2=dico_supers_hap2[dico_supers_hap2[0]==curated2].iloc[0,1]
        Paired.loc[index] = [super1, super2]
except Exception as e:
    print(f"An error occurred: {e}")

filtered_mashmap=pd.DataFrame(columns=mashmap_out.columns)

for index, row in mashmap_out.iterrows():
    if  ((Paired['Hap_1'] == row[[0,5]].iloc[0]) & (Paired["Hap_2"] == row[[0,5]].iloc[1])).any():
        filtered_mashmap.loc[len(filtered_mashmap)] = row

filtered_mashmap.to_csv(args.out_file, index=False, header=False, sep="\t")
