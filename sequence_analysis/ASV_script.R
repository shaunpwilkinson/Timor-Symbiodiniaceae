## DADA2 and insect classification pipeline script for Symbiodiniaceae NGS sequences

if(!require("insect")) install.packages("insect")
library(insect)
## set working directory
## sequences should be stored in subdirectory named 'seqs'
## create new directory 'dada2' containing new directory 'fastq'
## raw data available from open science framework (OSF) identifier ###
fpath <- "seqs/FASTQ_Generation_2017-11-28_18_33_52Z-64084105/"
samps <- dir(fpath)
## Took ~ 24h to write out FASTQ files with primers trimmed
for(i in seq_along(samps)){
  fpath2 <- paste0(fpath, samps[i])
  R1 <- readFASTQ(paste0(fpath2, "/", dir(fpath2)[1]))
  R2 <- readFASTQ(paste0(fpath2, "/", dir(fpath2)[2]))
  R1t <- trim(R1, up = "GTGAATTGCAGAACTCCGTG")
  R2t <- trim(R2, up = "CCTCCGCTTACTTATATGCTT")
  dif = " [12]"
  matches <- match(gsub(" [12]", "", names(R1t)), gsub(dif, "", names(R2t)))
  discards <- is.na(matches) 
  R1out <- R1t[!discards]
  matches <- matches[!discards]
  R2out <- R2t[matches]
  writeFASTQ(R1out, file = paste0("dada2/fastq/", gsub("\\.gz$", "", dir(fpath2)[1])))
  writeFASTQ(R2out, file = paste0("dada2/fastq/", gsub("\\.gz$", "", dir(fpath2)[2])))
}
## install and load dada2 package
source("https://bioconductor.org/biocLite.R")
biocLite("dada2")
library("dada2")
setwd("dada2/")
path <- "fastq"
list.files(path)
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1:2) 
## original dada script produced duplicates
sample.names <- apply(sample.names, 2, paste0, collapse = "_")
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
filtFs <- file.path("filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("filtered", paste0(sample.names, "_R_filt.fastq.gz"))
# On Windows set multithread=FALSE for filterAndTrim
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(220,200),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) 
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
saveRDS(derepFs, file = "derepFs.rds")
saveRDS(derepRs, file = "derepRs.rds")
saveRDS(errF, file = "errF.rds")
saveRDS(errR, file = "errR.rds")
saveRDS(out, file = "out.rds")
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
saveRDS(dadaFs, file = "dadaFs.rds")
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
saveRDS(dadaRs, file = "dadaRs.rds")
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
## save work here
saveRDS(mergers, file = "mergers.rds")
seqtab <- makeSequenceTable(mergers)
dim(seqtab) # 381 3759
table(nchar(getSequences(seqtab)))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
                                    multithread=TRUE, verbose=TRUE)
# Identified 1127 bimeras out of 3759 input sequences.
dim(seqtab.nochim) #381 2632
sum(seqtab.nochim)/sum(seqtab) # ~ 1% chimeras overall
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), 
               sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: 
# e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
## save work here
saveRDS(seqtab.nochim, file = "seqtabnochim.rds")


## Next part done on Ubuntu 16.04 server with VSEARCH installed
library(insect)
## set working directory to folder containing dada outout
stab <- readRDS("seqtabnochim.rds")
seqs <- colnames(stab)
stab <- t(stab)
seqdna <- char2dna(seqs) ## converts to DNAbin object
## download taxonomic classifier and read into R
URL <- "https://github.com/shaunpwilkinson/symulator/raw/master/data/classification_tree2.rds"
download.file(URL, destfile = "classification_tree2.rds", mode = "wb")
tree <- readRDS("classification_tree2.rds")
## score sequences against Symbiodiniaceae profile HMM
model <- decodePHMM(attr(tree, "model"))
score <- function(x, m) aphid::Viterbi(m, x)$score
cl <- parallel::makeCluster(6)
scores <- parallel::parSapply(cl, seqdna, score, model)
parallel::stopCluster(cl)
## remove low-scoring sequences (non-symbionts)
discards <- scores < 50
stab <- stab[!discards, ]
seqdna <- seqdna[!discards] # 2522 ASVs remaining
names(seqdna) <- hash(seqdna)

## OTU clustering with VSEARCH
tf1 <- tempfile(fileext = ".fa")
tf2 <- tempfile(fileext = ".aln")
writeFASTA(seqdna, tf1)
system(paste0("vsearch --cluster_smallmem ", tf1, " --alnout ", tf2, " --id 0.97 --usersort"))
otus <- scan(tf2, what = "", sep = "\n", quiet = TRUE)
queries <- otus[grepl("^Query", otus)]
queries <- gsub("^Query >", "", queries)
targets <- otus[grepl("^Target", otus)]
targets <- gsub("^.+>", "", targets)
centers <- unique(targets)
queries <- c(queries, centers)
targets <- c(targets, centers)
singles <- names(seqdna)[!names(seqdna) %in% queries]
queries <- c(queries, singles)
targets <- c(targets, singles)
clusters <- insect:::.point(targets)
names(clusters) <- queries
clusters <- clusters[match(names(seqdna), names(clusters))]
centrals <- names(clusters) %in% unique(targets)

## trim sequences to align globally with model
seqdna <- shave(seqdna, model, cores = 6)
ohashes <- hash(seqdna)
# get rereplication indices
opointers <- insect:::.point(ohashes) 
seqdna <- seqdna[!duplicated(opointers)]
## classify sequences to get taxon IDs
classifs <- classify(seqdna, tree, cores = 6)
clades <- classifs$taxon[opointers] # rereplicate clades
clades[clades == "root"] <- "UK" # unknown clade
oscores <- classifs$score[opointers] # rereplicate scores
## output table
front <- data.frame(variant = seq_along(clades),
                    clade = clades,
                    score = oscores,
                    otu = clusters,
                    central = centrals,
                    stringsAsFactors = FALSE)
back <- data.frame(hash = paste(insect::hash(char2dna(rownames(stab)))),
                   sequence = rownames(stab),
                   stringsAsFactors = FALSE)
rownames(stab) <- NULL
out <- cbind(front, stab, back)
write.csv(out, file = "ASVs.csv", row.names = FALSE)
## done

