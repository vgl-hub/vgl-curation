#!/usr/bin/env python3

import function
import textwrap
import os
import argparse



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


    function.AGPcorrect(args.fasta, args.agp, args.output_dir+'corrected_agp')

    function.hap_split(args.output_dir+'corrected_agp',args.output_dir+'Hap_1/hap1.agp',args.output_dir+'Hap_2/hap2.agp')

    print("\nRun unloc on Haplotype 1\n")
    function.unloc(args.output_dir+'Hap_1/hap1.agp',args.output_dir+'Hap_1/')
    print("\nRun unloc on Haplotype 2\n")
    function.unloc(args.output_dir+'Hap_2/hap2.agp',args.output_dir+'Hap_2/')


if __name__ == "__main__":
    main()
        