library(ggplot2)
library(tidyverse)
library(stringr)

rm(list = ls())
# Load readr
library(readr)
library(dplyr)

#set to your folder that has overview.txt and dbsub.out files
setwd("/Users/tejasn/Documents/OMalley_lab/Code/output_fungi/Neocon1//")
organism = substr(getwd(),55,65)
# Read text file with headers
overview_df= read_tsv('overview.txt')



#rename a few columns
colnames(overview_df)[2]= "EC"
colnames(overview_df)[6]= "Number_of_Tools"
colnames(overview_df)[1]= "GeneID"

setwd("/Users/tejasn/Documents/OMalley_lab/Code/output_fungi/Public_Neo_CAZy_content/")

#comparison with mycocosm
mycocosm_df = read_tsv(file = 'Neosp1_FilteredModels5_annotations_cazy_2023-08-16.tab')


overviewsubset<- subset(overview_df, Number_of_Tools!="1")


#reorder overview file by GeneID
overviewsubset = overviewsubset[order(overviewsubset$GeneID),] 

#read in dbcan_sub file
dbcan_sub_df = read_tsv('dbsub.out')  
colnames(dbcan_sub_df)[6]= "GeneID"

#reorder dbcan_sub file
dbcan_sub_df =dbcan_sub_df[order(dbcan_sub_df$GeneID),] 

#merge, rename '-' to unassigned in substrate
overview_dbcansub_tools2 = merge(overviewsubset,dbcan_sub_df,"GeneID")
overview_dbcansub_tools2$Substrate[overview_dbcansub_tools2$Substrate=='-'] = 'unassigned'

#identify columns that are duplicated by geneID (column 1)
merged_duplicates_boolean = duplicated(overview_dbcansub_tools2[,1])

#list of duplicated rows
merged_duplicates = overview_dbcansub_tools2[merged_duplicates_boolean,]
#doublechecking to make sure duplicates satisfy number of tools>1
#merged_duplicates = merged_duplicates[order(merged_duplicates$Number_of_Tools),] 

#anti_join to have unique GeneID dbcan_sub_outputs_only
unique_GeneIDs = anti_join(overview_dbcansub_tools2,merged_duplicates)
colnames(unique_GeneIDs)[1]= "GeneID"
write_csv(unique_GeneIDs,'unique_protein_IDs_dbcan_two_or_more_tools.csv')


#counting overall length of CAZymes

hmmer_df = read_tsv('hmmer.out')
colnames(hmmer_df)[4]= "gene_length"
colnames(hmmer_df)[3]= "GeneID"
hmmer_genes_duplicates_boolean = duplicated(hmmer_df[,3])
hmmer_duplicates = hmmer_df[hmmer_genes_duplicates_boolean,]
hmmer_df_uniques = anti_join(hmmer_df,hmmer_duplicates)
hmm_df_subset= merge(hmmer_df_uniques,overviewsubset,"GeneID")
total_cazyme_length = sum(hmm_df_subset$gene_length)

#summing CAZyme types for all geneIDs with two or more tools

AA = sum(str_count(overviewsubset$HMMER, "AA"))

CE = sum(str_count(overviewsubset$HMMER, "CE"))

GH = sum(str_count(overviewsubset$HMMER, "GH"))

GT = sum(str_count(overviewsubset$HMMER, "GT"))

PL = sum(str_count(overviewsubset$HMMER, "PL"))


#cazyme types
CAZyme_vector = c(AA,CE, GH, GT, PL)
cazyme_sum = sum(CAZyme_vector)
CAZyme_vector = c(CAZyme_vector,cazyme_sum,total_cazyme_length)
CAZyme_type = c('AA','CE','GH','GT','PL','Sum','total_cazyme_length bp')
CAZyme_df = data.frame(CAZyme_type,CAZyme_vector)
colnames(CAZyme_df)= c('CAZyme_type','genome_count')
cazyme_type_filename = str_c('/Users/tejasn/Documents/OMalley_lab/Code/output_fungi/',organism,'cazymetypetable.csv')
write_excel_csv(CAZyme_df,cazyme_type_filename)


fraction_CAZyme = CAZyme_df$genome_count/cazyme_sum
CAZyme_df$fraction = fraction_CAZyme

#graphing CAZymes
ggplot(CAZyme_df, aes(x = CAZyme_type,y=genome_count))+
  geom_bar(stat='identity')+  coord_flip() + ggtitle(organism)

#fraction CAZymes
#ggplot(CAZyme_df, aes(x = CAZyme_type,y=fraction))+
 # geom_bar(stat='identity')+  coord_flip() + ggtitle(organism)


#use overview_dbcansub_tools2 dataframe because dbcan-sub output 
#has one CAZyme per row (and predicted substrates),
#and some geneIDs have multiple CAZymes

#types of substrates:
#use R to make list of substrates from output columns


#getting rid of commas and making a separate list
substrates = gsub(',',' ', overview_dbcansub_tools2$Substrate)
head(substrates)
write.csv(substrates,file = 'testlist3.csv')

sum(substrate_bool=unlist(strsplit(substrates, " "))!='' ,na.rm=TRUE)

testlist = unlist(strsplit(substrates, " "))
write.csv(testlist,file = 'testlist.csv')



#getting all uniques and taking out blanks and -, optional
substrate_list = unique(unlist(strsplit(substrates, " ")))
substrate_list = substrate_list[substrate_list!='']
#substrate_list_without_unassigned = substrate_list[substrate_list!='unassigned']
length(substrate_list)
#NOTE THAT SUBSTRATE_LIST includes duplicated values for one enzyme.
#grep fixes this issue in for loop a few lines down


#making a dataframe with substrate_count and initializing with zeroes
substrate_count = data.frame(substrate_list,integer(length(substrate_list)))
colnames(substrate_count)=c('substrate','occurrences')



#initialize for loop that counts for each name in the substrate list
i=1
substrate_list[2]
for (substrate in substrate_list)
{
  print(substrate)
  substrate_count$occurrences[i] = length(grep(substrate_list[i],overview_dbcansub_tools2$Substrate))
   i=i+1
  
}

#computing fractional occurrence
sum_substrate_count = sum(substrate_count$occurrences)
sum_substrate_count
fraction_substrate = substrate_count$occurrences/sum_substrate_count
substrate_count$fraction = fraction_substrate

substrate_filename = str_c('/Users/tejasn/Documents/OMalley_lab/Code/output_fungi/',organism,'_substrate_table.csv')

write_excel_csv(substrate_count,substrate_filename)
#graphing occurrences
ggplot(substrate_count, aes(x = substrate,y=occurrences))+
  geom_bar(stat='identity')+  coord_flip() + ggtitle(organism)

#plotting without unassigned

ggplot(subset(substrate_count,subset=substrate_count$substrate!='unassigned')
, aes(x = substrate,y=occurrences))+
  geom_bar(stat='identity')+  coord_flip() + ggtitle(organism)


#graphing fraction
#ggplot(substrate_count, aes(x = substrate,y=fraction))+
 # geom_bar(stat='identity')+  coord_flip() + ggtitle(organism)
