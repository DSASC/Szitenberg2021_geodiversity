
library(Matrix);
library(dada2);

path <- "~/GeoDiversity/temp" # CHANGE ME to the directory containing the fastq files after unzipping.
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
# Place filtered files in filtered/ subdirectory
path1 = "~/GeoDiversity/cut_reads"                                                        
filtFs <- file.path(path1, "ITS1", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path1, "ITS1", paste0(sample.names, "_R_filt.fastq.gz"))                                                        
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(125,125),
              maxN=0, maxEE=c(4,5), truncQ=10, rm.phix=TRUE, trimLeft=0,
              compress=TRUE, multithread=TRUE)
# learn F errors
errF <- learnErrors(filtFs, multithread=TRUE)
print('learn R errors')
errR <- learnErrors(filtRs, multithread=TRUE)

# derep F
derepFs <- derepFastq(filtFs, verbose=TRUE)
# derep R
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# dada F
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
# dada R
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# merging
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, justConcatenate=TRUE)
seqtab <- makeSequenceTable(mergers)

# chimera removal                                                    
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# assign taxonomy                                                        
taxa <- assignTaxonomy(seqtab, "~/GeoDiversity/DB/sh_general_release_dynamic_02.02.2019.fasta")

# write outputs                                                   
uniquesToFasta(seqtab, 'dada/ITS1/rep-seqs.fasta')
write.csv(seqtab, file = "dada/ITS1/table.csv")
write.csv(taxa, file = "dada/ITS1/taxa.csv")
