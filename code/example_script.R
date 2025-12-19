library(devtools)
library(stringr)
source('/Users/zuckem04/Documents/Side_Projects/Folate/pipeline_scripts/signature_functions.R')

#example vcfs:
vcfs <- list.files('/Users/zuckem04/Documents/Side_Projects/Folate/strelka_vcfs_new', full.names = T)
sampleNames <- str_split_fixed()

#Load cosmic mutation dictionary:
cosmic <- read.table('/Users/zuckem04/Documents/Side_Projects/Folate/COSMIC_v3.4_SBS_GRCh37.txt', sep='\t', header=T)

#Getting signatures:
signatures <- assignSignatures(vcfs, dictionary=cosmic)
exposures <- signatures$contribution
mut_mat <- signatures$mut_table

#plotting/making mut profile
mut_profile <- mutationProfilePlot(mut_mat, sampleName=colnames(mut_mat)[2])

#plotting exposures:
mutSigPlot2(exposures)


