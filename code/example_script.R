library(devtools)
library(stringr)
source('/Users/zuckem04/Documents/Side_Projects/Folate/pipeline_scripts/signature_functions.R')

#example vcfs:
vcfs <- list.files('/Users/zuckem04/Documents/Side_Projects/Folate/strelka_vcfs_new', full.names = T)
sampleNames <- str_split_fixed(vcfs, "/", Inf)[,8]

#Load cosmic mutation dictionary:
cosmic <- read.table('/Users/zuckem04/Documents/Side_Projects/Folate/COSMIC_v3.4_SBS_GRCh37.txt', sep='\t', header=T)

#Getting signatures:
signatures <- assignSignatures(vcfs, dictionary=cosmic, sampleNames = sampleNames)
exposures <- signatures$contribution
mut_mat <- as.data.frame(signatures$mut_table)
colnames(mut_mat) <- sampleNames
write.csv(mut_mat, file="/Users/zuckem04/Documents/temp/merged_mutation_profiles.csv")

#plotting/making mut profile
mut_profile <- mutationProfilePlot(mut_mat, sampleName=sampleNames[1], 
                                   savePath='/Users/zuckem04/Documents/temp/sig_examples')

#run on list of samples
for(i in 1:length(sampleNames)){
  mutationProfilePlot(mut_mat, sampleName=sampleNames[i], 
                                     savePath='/Users/zuckem04/Documents/temp/sig_examples')
}

#plotting exposures:
mutSigPlot2(exposures)


