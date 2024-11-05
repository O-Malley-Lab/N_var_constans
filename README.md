# N_var_constans
genomic, transcriptomic, functional analysis of neocallimastix var constans

genome-level files: 
Neocon1_FilteredModels1_annotations_cazy_2023-08-16.tab contains Mycocosm-filtered and annotated CAZymes in N. sp. Constans.

Neocon1_GeneCatalog_proteins_20220615.aa.fasta is an amino acid FASTA file of filtered, translated ORFs in N. sp. Constans, also available on Mycocosm.org

PF02013.hmm is the hmmer profile for dockerin, from PFAM

cohesin3.hmm.txt is the hmmer profile for cohesin, which are generally present in multiple copies in fungal scaffoldins. This file is an updated version of the cohesin hmm file published with http://dx.doi.org/10.1038/nmicrobiol.2017.87 

neocon1_run_dbcan_parser.R is an R script to parse dbcan3 output

unique_protein_IDs_dbcan_two_or_more_tools.csv is the filtered list of CAZymes from dbcan3 output

