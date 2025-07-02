import pandas as pd 
import argparse
from io import StringIO
from natsort import natsorted


parser = argparse.ArgumentParser(
                    prog='filter_mashmap_with_tagged_pairs',
                    description='Filter Mashmap output to keep only Scaffolds paired with the tags Micro1 and Micro2',
                    usage='filter_mashmap_with_tagged_pairs.py -1 Hap_1/ -2 Hap_2/ -m mashmap.out -o filtered_mashmap.out')
 
parser.add_argument('-1', '--hap1', dest="hap1",required=True, help='Path to the chromosome assignment file for  Haplotype 1 (inter_chr.tsv) ')  
parser.add_argument('-2', '--hap2', dest="hap2", required=True, help='Path to the chromosome assignment file for Haplotype 2 (inter_chr.tsv) ') 
parser.add_argument('-q', '--query', dest="query", default="Hap_2", help='Haplotype use as Query for MashMap: Hap_1 or Hap_2 (Default Hap_2)')  
parser.add_argument('-r', '--reference', dest="reference", default="Hap_1", help='Haplotype use as reference for MashMap: Hap_1 or Hap_2 (Default Hap_1)')  
parser.add_argument('-a', '--agp', dest="agp", required=True, help='Path to the curated AGP file (containing the tags Micro1 and Micro2)')       
parser.add_argument('-m', '--mashmap', dest="mashmap", required=True, help='Path to the curated Mashmap output')  
parser.add_argument('-o', '--out', dest="out_pre", required=True, help='Path to output Prefix')  

args = parser.parse_args()


if args.query!="Hap_1" and args.query!="Hap_2":
    raise ValueError("The query parameter should be Hap_1 or Hap_2")
if args.reference!="Hap_1" and args.reference!="Hap_2":
    raise ValueError("The reference parameter should be Hap_1 or Hap_2")

try:
    dico_supers_hap1= pd.read_csv(args.hap1, sep="\t", header=None)
    dico_supers_hap2= pd.read_csv(args.hap2, sep="\t", header=None)
except FileNotFoundError:
    print(f"Error: The file inter_chr.tsv was not found.")

try:
    mashmap_out=pd.read_csv(args.mashmap, sep="\t", header=None)
except FileNotFoundError:
    print(f"Error: The file '{args.mashmap}' was not found.")

mashmap_out=mashmap_out.rename(columns={0:args.query, 1:'Query length',2:'Query 0-based start',3:'Query 0-based end',4:'Strand',5:args.reference, 6:'Target length', 7:'Target 0-based start', 8:'Target 0-based end',12:'Mapping nucleotide identity'})

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

filtered_agp = pd.read_csv(StringIO(output_agp),sep='\t',header=None)

dico_micro=pd.DataFrame(columns=['Hap_1', 'Hap_2'])

even_value = None
even_hap=None

for index, row in filtered_agp.iterrows():  
    if index % 2 == 0:
        even_value=row.iloc[0] 
        even_hap=row.iloc[10]
    else:
        index_new=len(dico_micro)
        tmp=pd.DataFrame(columns=['Hap_1', 'Hap_2'])
        tmp.loc[0,row.iloc[10]]=row.iloc[0]
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
            raise ValueError("Names from the curated agp were not founds in the hap.unlocs.no_hapdups.agp file. Verify the names are consistants.")
        super1=dico_supers_hap1[dico_supers_hap1[0]==Scaffold1].iloc[0,1]
        super2=dico_supers_hap2[dico_supers_hap2[0]==Scaffold2].iloc[0,1]
        Paired.loc[index] = [super1, super2]
except Exception as e:
    print(f"An error occurred: {e}")

filtered_mashmap=pd.DataFrame(columns=mashmap_out.columns)

for index, row in mashmap_out.iterrows():
    if  ((Paired['Hap_1'] == row[[args.query,args.reference]].iloc[1]) & (Paired["Hap_2"] == row[[args.query,args.reference]].iloc[0])).any():
        filtered_mashmap.loc[len(filtered_mashmap)] = row

orientations = filtered_mashmap[[args.query,"Strand",args.reference]]

# Group by Reference and Query, then count Sign occurrences
counts = orientations.groupby(["Hap_1","Hap_2"])["Strand"].value_counts()

# Get the sign with the highest count per group
most_frequent_sign = counts.groupby(level=[0, 1]).idxmax().apply(lambda x: x[2])  # x is a tuple (Reference, query, sign)

# Combine into a final dataframe
result = most_frequent_sign.to_frame(name='Main Orientation')

result=result.reset_index()
result['Hap_1'] = pd.Categorical(result['Hap_1'], categories=natsorted(result['Hap_1'].unique()), ordered=True)
result = result.sort_values(by='Hap_1')

result.to_csv(args.out_pre+"orientation.tsv", index=True, header=True, sep="\t")

to_reverse = pd.DataFrame(columns=['Command', 'Hap2_SUPER'])
to_reverse['Hap2_SUPER']=result[result['Main Orientation']=='-']["Hap_2"]
to_reverse['Command']="RVCP"


to_reverse.to_csv(args.out_pre+"rvcp.sak", index=False, header=False, sep="\t")

renaming=result[['Hap_1','Hap_2']]
renaming=renaming.rename(columns={'Hap_1':'New_name','Hap_2':'Old_name'})
renaming['action']="RENAME"
renaming[['action','Old_name','New_name']].to_csv(args.out_pre+"renaming_hap2.tsv", index=False, header=True, sep="\t")

