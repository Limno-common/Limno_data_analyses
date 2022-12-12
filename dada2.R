library("stringr")
library(dada2) 
packageVersion("dada2")

setwd("/path/to/working/directory/with/fastq/files")
list.files() 
forward_reads <- sort(list.files(pattern="_R1_001.fastq.gz")) 
reverse_reads <- sort(list.files(pattern="_R2_001.fastq.gz")) 

#Visually check the quality of your reads. Is there any big drop in quality at the beginning or the end? A gradual drop in quality at the end ('3) is normal
#one plot for each file, you can look at as many as you want
plotQualityProfile(forward_reads[1:5])
plotQualityProfile(reverse_reads[1:5]) 

#Create more usable names for your samples
sample.names <- sapply(strsplit(basename(forward_reads), "_"), `[`, 1)
sample.names <- str_remove(sample.names, "UC-2885-") #change according to what part of the file name you want to remove
sample.names

filtFs <- file.path("/path/to/working/directory", "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("/path/to/working/directory", "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#Quality filtering of your reads
#dada2 does not allow any Ns in sequence, you therefore have to remove these bits, otherwise you'll loose the entire read
#trimLeft: remove primers, potential linking sequences and bad quality bases at the 5'end from fwd and rev reads, you can trim different amounts of bases for fwd and rev
#truncQ: Q2 = # in Illumina -> 63% chance of erroneous base assignment, bases until this quality threshold will be removed from 3'
#maxEE: maximum number of “expected errors” allowed in a read, can be different for fwd and rev
out <- filterAndTrim(forward_reads, filtFs, reverse_reads, filtRs,
                     rm.phix=TRUE, 
                     multithread=TRUE,
                     compress=TRUE,
                     trimLeft = c(21,54), 
                     truncQ = 2,
                     maxEE=c(3,6), 
                     minLen = 235)
#always good to save your R objects in between
saveRDS(out, file = "FilterOut.rds")
#out <- readRDS("FilterOut.rds")
head(out)

#Learn the error rate in your reads to correct for it later when calling ASVs
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#ASV calling
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
saveRDS(dadaFs, file = "dadaFs.rds")
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
saveRDS(dadaRs, file = "dadaRs.rds")

#dadaRs <- readRDS("dadaRs.rds")
#dadaFs <- readRDS("dadaFs.rds")

#merge fwd and rev reads
#default minOverlap is 20bp, choose something appropriate for your read length and expected sequence length, default maxMismatch is 0
merge <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap = 5, maxMismatch = 1, verbose=TRUE) head(merge[[1]])

seqtab <- makeSequenceTable(merge)
dim(seqtab)

#Check the length of all merged sequences. What length did you expect? Which sequences are too short to be true?
table(nchar(getSequences(seqtab)))
#decide which sequences to keep
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 480:522]
table(nchar(getSequences(seqtab2)))

#Remove chimera
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)

dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab2)

#Check how many reads were lost at each step
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(merge, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#Taxonomic assignment of ASVs based on some reference database. I use PR2 (https://github.com/pr2database/pr2database/releases), but you can also use SILVA, greengenes, etc.
taxa <- assignTaxonomy(seqtab.nochim, "path/to/reference/database/file", multithread=TRUE, taxLevels = c("Kingdom","Supergroup", "Phylum", "Class","Order","Family","Genus","Species"))
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#giving our seq headers more manageable names (ASV_1, ASV_2...) 
asv_seqs <- colnames(seqtab.nochim) 
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character") 
for (i in 1:dim(seqtab.nochim)[2]) { 
  asv_headers[i] <- paste(">ASV", i, sep="_") 
} 
asv_headers

#create fasta file with sequence of each ASV
asv_fasta <- c(rbind(asv_headers, asv_seqs)) 
write(asv_fasta, "path/to/folder/where/you/want/to/save/the/file")

#create ASV count table 
asv_tab <- t(seqtab) 
row.names(asv_tab) <- sub(">", "", asv_headers) 
write.table(asv_tab, "path/to/folder/where/you/want/to/save/the/file", sep="\t", quote=F) 

#create taxonomic table 
asv_tax <- taxa 
row.names(asv_tax) <- sub(">", "", asv_headers) 
write.table(asv_tax, "path/to/folder/where/you/want/to/save/the/file", sep="\t", quote=F)


