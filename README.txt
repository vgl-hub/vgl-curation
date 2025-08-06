================================
Process Curated Genome Assembly
================================

ProcessCurated is a set of tools to process curated genome assemblies. In combination with gfastats and mashmap, tt reconciliate the AGP file manually curated in PretextView to rename, reorient, and sort the assemblies to get them ready for submission. 


Available tools:

* AGPCorrect : Correcting AGP for sequence lengths

* hap_split: Splitting the haplotypes from the corrected AGP

* unloc: Assigning unlocs before the agp is imposed on the fasta

* chromosome_assignment: Substituting scaffold for chromosome assignments

* sak_generation: Processing Mashmap output for reorientation and renaming of scaffolds in haplotype 2.


Tools usage
===========

AGPCorrect
-------------

Inputs:

* Fasta file of the Assembly with the two haplotypes

* Curated AGP generated with PretextView (See manual in `docs` for curation tips)


Usage::

    AGPcorrect assembly.fasta curated.agp > corrected.agp 

Output:

* Corrected AGP file


hap_split
-------------

Inputs:

* Corrected AGP file (-a)

* Output AGP file for Haplotype 1 (-1) 

* Output AGP file for Haplotype 2 (-2) 


Usage::

    hap_split -a corrected.agp -1 Hap_1/hap1.agp -2 Hap_2/hap2.agp 


Output:

* AGP file for each haplotype


Run whole analysis
==================




Contributors
============