library(MutationalPatterns)
library(ref_genome, character.only = TRUE)

assignSignatures <- function(vcfs, sampleNames=NULL, ref="BSgenome.Hsapiens.UCSC.hg19", dictionary="cosmic_snv"){
  ref_genome <- ref
  if(is.null(sampleNames)){samNanes <- vcfs}else{samNames <- sampleNames}
  grl <- read_vcfs_as_granges(vcfs_finished, samNames, ref_genome, type="snv")
  mut_mat <- mut_matrix(vcf_list = grl_finished, ref_genome = ref_genome)
  fit_res <- fit_to_signatures(mut_mat, cosmic)
}

#Load cosmic mutation dictionary:
cosmic <- read.table('~/projects/folate/cosmic/COSMIC_v3.4_SBS_GRCh37.txt', sep='\t', header=T)
rownames(cosmic) <- cosmic$Type
cosmic <- as.matrix(cosmic[,2:ncol(cosmic)])
cosmic <- cosmic[match(rownames(mut_mat_finished), rownames(cosmic)),]

#de novo extraction:
#...



