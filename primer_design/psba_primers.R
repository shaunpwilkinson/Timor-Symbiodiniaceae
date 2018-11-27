#####################################################################################
############### PRIMER DESIGN SCRIPT FOR UNIVERSAL PSBANCR PRIMERS ##################
#####################################################################################

library(insect)

## download reference sequences from GenBank
rdb <- searchGB(paste0("Symbiodiniaceae[ORGN]+AND+PSBA[GENE]+NOT+predicted[TITL]",
                       "+NOT+unverified[TITL]+NOT+environmental[ORGN]+AND+100:5000[SLEN]"))
taxids <- as.integer(sub(".+\\|", "", names(rdb)))
tdb <- taxonomy() ## downloaded Nov 2018
tdb <- prune_taxonomy(tdb, taxids)
saveRDS(rdb, file = "rdb.rds")
saveRDS(tdb, file = "primer_design/tdb.rds")

## read in databases
rdb <- readRDS("primer_design/rdb.rds") # reference database
tdb <- readRDS("primer_design/rdb.rds") # taxonomy map
SYMPSBANCRF <- "AATCGTGCTGATCTAGGWATGG"
SYMPSBANCRR <- "GAGACGATTTGTTGTGGATAG"

## trim LHS of primer bind region
forw_seqs <- virtualPCR(rdb, up = SYMPSBANCRF, trimprimers = FALSE, cores = 6) #44
rev_seqs <- virtualPCR(rdb, up = rc(SYMPSBANCRR), trimprimers = FALSE, cores = 6) #77

## only keep sequences for which both primer bind regions are available
keeps <- names(forw_ali) %in% names(rev_ali) #19
forw_ali <- forw_seqs[keeps]
keeps <- names(rev_ali) %in% names(forw_ali)
rev_ali <- rev_seqs[keeps]
rev_ali <- rev_ali[names(forw_ali)]

## trim RHS
forw_ali <- shave(forw_ali, char2dna(SYMPSBANCRF))
rev_ali <- shave(rev_ali, char2dna(rc(SYMPSBANCRR)))

## remove duplicates
concats <- paste0(sub(".+\\|", "", names(forw_ali)), dna2char(forw_ali), dna2char(rev_ali))
discards <- duplicated(concats) #10
forw_ali <- forw_ali[!discards]
rev_ali <- rev_ali[!discards]

## build alignments
forw_ali <- aphid::align(forw_ali)
rev_ali <- aphid::align(rev_ali)

## find lineage metadata
taxids <- as.integer(sub(".+\\|", "", rownames(forw_ali)))
lineages <- get_lineage(taxids, tdb)
unname(sapply(lineages, tail, 1))
sub("\\|.+", "", rownames(forw_ali))

## print alignments
matrix(apply(forw_ali, 1, dna2char), ncol = 1)
matrix(apply(rev_ali, 1, dna2char), ncol = 1)

## find positions in reference sequence
HG515017_CDS <- structure(rdb["HG515017|431157"][[1]][1:1029], class = "DNAbin")
match(1, aphid::Viterbi(HG515017_CDS, char2dna(SYMPSBANCRF), type = "semiglobal")$path) #958
match(1, aphid::Viterbi(HG515017_CDS, char2dna(rc(SYMPSBANCRR)), type = "semiglobal")$path) #58

## add missing clades D, E, F manually if possible
ftaxids <- as.integer(sub(".+\\|", "", names(forw_seqs)))
flins <- get_lineage(ftaxids, tdb)
unname(sapply(flins, tail, 1))
rtaxids <- as.integer(sub(".+\\|", "", names(rev_seqs)))
rlins <- get_lineage(rtaxids, tdb)
unname(sapply(rlins, tail, 1))

#####################################################################################