library(data.table)
library(ggplot2)
library(stringr)
library(ggpubr)
library(vcfR)
library(MutationalPatterns)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

#Some plotting functions:
mutSigPlot <- function(exposures, sample_order = NULL, main=NULL, showSignatures=NULL, topN=NULL, colorset=NULL, exposure=TRUE,
                       ylabl=NULL, xaxisText=TRUE, ylm=NULL){
  #Show top n signatures:
  #Note: forcing topN to be 73 since that's the maximum number of colors I can show:
  if(!is.null(topN)){sigs <- rownames(exposures[order(rowSums(exposures), decreasing = T),])[1:topN]}else{
    sigs <- rownames(exposures[order(rowSums(exposures), decreasing = T),])[1:73]}
  
  #if showSignatures isn't NULL, everything not in these becomes 'other:
  if(!is.null(showSignatures)){
    sigs <- showSignatures
    expo <- rbind(exposures[which(rownames(exposures) %in% sigs),], 
                  colSums(exposures[which(!rownames(exposures) %in% sigs),]))
    rownames(expo)[nrow(expo)] <- 'Other'
  }else{
    expo <- exposures}
  
  #colors:
  if(is.null(colorset)){
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    colorvec <- sample(col_vector, nrow(expo))
    names(colorvec) <- rownames(expo)
    colorvec[which(colorvec == 'Other')] <- 'gray'
  }else{colorvec <- c(colorset,'Other'='gray')}
  
  df <- Reduce(rbind, lapply(1:ncol(expo),function(i){
    data.frame('Signature'=rownames(expo),'Sample'=colnames(expo)[i],'Exposure'=expo[,i])}))
  if(!is.null(sample_order)){df$Sample <- factor(df$Sample, levels = sample_order)}
  if(is.null(main)){titl <- ''}else{titl <- main}
  
  #Restricting color vector to signatures present:
  colorvec <- colorvec[which(names(colorvec) %in% df$Signature)]
  if('unassigned' %in% df$Signature){colorvec <- c(colorvec,'unassigned'='black')}
  #Ordering:
  df$Signature <- factor(df$Signature, levels = unique(df$Signature))
  
  #Defining muts per sample:
  df$TMB <- sapply(1:nrow(df), function(i){sum(df$Exposure[which(df$Sample == df$Sample[i])])})
  
  #Some plotting options:
  if(xaxisText){xax <- element_text(angle = 90, vjust = .5, hjust=1)}else{xax <- element_blank()}
  if(!is.null(ylabl)){}else{if(exposure){ylabl <- 'Exposure'}else{ylabl <- 'Fraction'}}
  if(is.null(main) | main==''){theme_main <- element_blank()}else{theme_main <- element_text(hjust=.5)}
  if(is.null(ylm)){ymax <- max(df$TMB)}else{ymax <- ylm}
  
  if(exposure == TRUE){
    ggplot(data=df, aes(x=Sample, y=Exposure, fill=Signature)) + geom_bar(stat = 'identity', position = 'stack') + 
      theme_bw() + theme(axis.text.x = xax, plot.title = theme_main) + ylim(c(0,ymax)) + 
      xlab('') + ylab(ylabl) + scale_fill_manual(values = colorvec) + ggtitle(main)
  }else{
    df$Fraction <- sapply(1:nrow(df), function(i){df$Exposure[i]/df$TMB[i]})
    ggplot(data=df, aes(x=Sample, y=Fraction, fill=Signature)) + geom_bar(stat = 'identity', position = 'stack') + 
      theme_bw() + theme(axis.text.x = xax, plot.title = theme_main) + 
      xlab('') + ylab(ylabl) + scale_fill_manual(values = colorvec) + ggtitle(main)
  }
}

mutSigPlot2 <- function(exposure){
  expo <- exposure
  expo <- cbind('signature'=row.names(expo), expo)
  melted <- reshape2::melt(as.data.frame(expo), id.vars='signature')
  colnames(melted) <- c('Signature','Sample','Exposure')
  colors <- c('SBS18'='blue','SBS17a'='darkgreen','SBS17b'='lightgreen','SBS4'='purple')
  melted$Exposure <- as.numeric(melted$Exposure)
  gg <- ggplot(data=melted, aes(x=Sample, y=Exposure, fill=Signature)) + geom_bar(stat = 'identity', position = 'stack') + 
    theme_bw() + theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=0)) + 
    xlab('') + ylab('Exposure') + ggtitle("") + scale_fill_manual(values = colors)
  print(gg)
}

mutationProfilePlot <- function(SNV_catalogues, sampleName, savePath=NULL, main=''){
  if('data.frame' %in% class(SNV_catalogues)){
    temp <- data.frame('muttype'=rownames(SNV_catalogues), 'frequency'=SNV_catalogues[,sampleName])
  }else{
    temp <- data.frame('muttype'=names(SNV_catalogues), 'frequency'=SNV_catalogues)
  }
  spl <- str_split_fixed(temp$muttype, "", Inf)
  temp$Group <- str_c(spl[,3], spl[,4], spl[,5])
  temp$Group <- factor(temp$Group, levels = unique(temp$Group))
  temp$muttype <- factor(temp$muttype, levels = unique(temp$muttype))
  muttypePlot <- ggplot(data=temp, aes(x=muttype, y=frequency, fill=Group)) + geom_bar(stat = 'identity') + theme_bw() +
    theme(axis.text.x = element_text(angle=90, vjust=.5)) + xlab('') + ylab('Mutation Freq.') + ggtitle(main) + 
    theme(plot.title = element_text(hjust=.5)) + 
    scale_fill_manual(values = c('C>A'='#55B9E8','C>G'='black','C>T'='#D03D32','T>A'='#CAC8C9','T>C'='#AACD70','T>G'='#E6C7C5'))  
  if(!is.null(savePath)){
    pdf(paste0(savePath,'/',sampleName,'.pdf'), width = 15, height = 5)
    print(muttypePlot)
    dev.off()
    write.csv(temp, file=paste0(savePath,'/',sampleName,'.csv'),col.names=T, row.names=F, quote=F)
  }else{
    print(muttypePlot)
  }
  temp
}

cosmicPlot <- function(cosmic, signature, anno=""){
  #Plotting:
  #cosmic (for reference):
  if('Type' %in% colnames(cosmic)){
    cosmic_sig <- as.data.frame(cbind('MutationType'=as.data.frame(cosmic)[,'Type'], 
                                      'frequency'=as.data.frame(cosmic)[,signature]))
  }else{
    cosmic_sig <- as.data.frame(cbind('MutationType'=rownames(cosmic), 
                                      'frequency'=as.data.frame(cosmic)[,signature]))
  }
  cosmic_sig$frequency <- as.numeric(cosmic_sig$frequency)
  spl.cosmic <- str_split_fixed(cosmic_sig$MutationType, "", Inf)
  cosmic_sig$Group <- str_c(spl.cosmic[,3], spl.cosmic[,4], spl.cosmic[,5])
  cosmic_sig$Group <- factor(cosmic_sig$Group, levels = unique(cosmic_sig$Group))
  cosmic_sig <- cosmic_sig[order(cosmic_sig$Group),]
  cosmic_sig$MutationType <- factor(cosmic_sig$MutationType, levels = unique(cosmic_sig$MutationType))
  cosmicPlot <- ggplot(data=cosmic_sig, aes(x=MutationType, y=frequency, fill=Group)) + geom_bar(stat = 'identity') + theme_bw() + 
    theme(axis.text.x = element_text(angle=90, vjust=.5,hjust=1)) + xlab('') + ylab('Mutation Prob.') + 
    scale_fill_manual(values = c('C>A'='#55B9E8','C>G'='black','C>T'='#D03D32','T>A'='#CAC8C9','T>C'='#AACD70','T>G'='#E6C7C5')) + 
    theme(plot.title = element_text(hjust=.5)) + ggtitle(paste0(signature,anno))
  cosmicPlot
}


#Generate variant tables:
#currently supports strelka and sage
makeVarTable <- function(vcf, caller='strelka'){
  fix <- vcf@fix
  temp <- as.data.frame(fix[,c(1,2,4,5)])
  spl.ref <- str_split_fixed(temp$REF, '', Inf)
  spl.alt <- str_split_fixed(temp$ALT, '', Inf)
  #Annotating whether each variant is an SNV:
  temp$SNV <- sapply(1:nrow(spl.ref),function(i){length(which(spl.ref[i,] != '')) == 1 & length(which(spl.alt[i,] != '')) == 1})
  snvs <- temp[which(temp$SNV),]
  
  #extracting allele counts from gt object:
  gt <- as.data.frame(vcf@gt)[which(temp$SNV),]
  format <- as.data.frame(str_split_fixed(gt$FORMAT, ':', Inf))
  if(caller == 'strelka'){
    if(is.null(gt$tumor)){
      tumor <- as.data.frame(str_split_fixed(gt[[length(gt)]], ':', Inf))
    }else{
      tumor <- as.data.frame(str_split_fixed(gt$TUMOR, ':', Inf)) 
    }
    colnames(format) <- colnames(tumor) <- format[1,]
    colnames(tumor)[5:8] <- c('A','C','G','T')
    #Filtered depth:
    snvs$counts_filtered <- as.numeric(tumor$DP) - as.numeric(tumor$FDP)
    snvs$ref_count <- as.numeric(sapply(1:nrow(snvs), function(i){strsplit(tumor[i,which(colnames(tumor) == snvs$REF[i])], ',')[[1]][1]}))
    snvs$alt_count <- as.numeric(sapply(1:nrow(snvs), function(i){strsplit(tumor[i,which(colnames(tumor) == snvs$ALT[i])], ',')[[1]][1]}))
    snvs$total_count <- rowSums(sapply(5:8, function(i){as.numeric(str_split_fixed(tumor[,i], ',', 2)[,1])}))
  }else{
    if(is.null(gt$tumor)){
      tumor <- as.data.frame(str_split_fixed(gt[[length(gt)]], ':', Inf))
    }else{
      tumor <- as.data.frame(str_split_fixed(gt$tumor, ':', Inf)) 
    }
    colnames(tumor) <- unlist(format[1,])
    #Filtered depth:
    counts <- str_split_fixed(tumor$AD, ',', 2)
    snvs$counts_filtered <- as.numeric(counts[,1]) + as.numeric(counts[,2])
    snvs$ref_count <- as.numeric(counts[,1])
    snvs$alt_count <- as.numeric(counts[,2])
    snvs$total_count <- NA
  }
  #Changing a couple column names:
  colnames(snvs)[1:2] <- c('chr','position')
  #Numericizing:
  for(j in c(2,6:ncol(snvs))){snvs[,j] <- as.numeric(as.character(snvs[,j]))}
  #adding vaf:
  snvs$vaf <- snvs$alt_count/(snvs$alt_count + snvs$ref_count)
  #Should do some sanity checks here
  snvs
}

#dictionary will be cosmic usually. 
assignSignatures <- function(vcfs, dictionary, sampleNames=NULL, ref="BSgenome.Hsapiens.UCSC.hg19"){
  ref_genome <- ref
  dict <- dictionary
  if('Type' %in% colnames(dict)){
    rownames(dict) <- dict$Types
    dict <- dict[,2:ncol(dict)]
  }
  #must be in matrix format
  dict <- as.matrix(dict)
  if(is.null(sampleNames)){samNames <- vcfs}else{samNames <- sampleNames}
  
  #read in vcfs as granges objects
  grl <- read_vcfs_as_granges(vcfs, samNames, ref_genome, type="snv")
  mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
  
  #run signatures
  fit_res <- fit_to_signatures(mut_mat, dict)
  fit_res$mut_table <- mut_mat 
  fit_res
}

#de novo extraction:
extract_snv_signatures <- function(vcfs, rank, sampleNames=NULL, ref="BSgenome.Hsapiens.UCSC.hg19"){
  if(is.null(sampleNames)){samNames <- vcfs}else{samNames <- sampleNames}
  grl <- read_vcfs_as_granges(vcfs, samNames, ref_genome, type="snv")
  mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
  mut_mat <- mut_mat + 0.0001
  nmf_res <- extract_signatures(mut_mat, rank = rank, nrun = 10, single_core = TRUE)
  nmf_res
}


