library(insect)

setwd("~/Dropbox/East_Timor/Rutherford_Project/Part3_treecongruence/Josh_Paper1/R/PSBA/")
rdb <- searchGB(paste0("Symbiodiniaceae[ORGN]+AND+PSBA[GENE]+NOT+predicted[TITL]",
                       "+NOT+unverified[TITL]+NOT+environmental[ORGN]+AND+100:5000[SLEN]"))
saveRDS(rdb, file = "rdb.rds")
## reference sequence - coding sequence postion looked up manually

rdb <- readRDS("rdb.rds")
tdb <- readRDS("~/Dropbox/R/insect/taxonomy/tdb_full.rds") # Nov 2018 version

psbA_1.0 <- "CWGTAGATATTGATGGWATAAGAGA"
psbA_3.0 <- "TTGAAAGCCATTGTYCTTACTCC"

SYMPSBANCRF <- "AATCGTGCTGATCTAGGWATGG"
SYMPSBANCRR <- "GAGACGATTTGTTGTGGATAG"

test <- virtualPCR(rdb, up = SYMPSBANCRF, down = SYMPSBANCRR, trimprimers = FALSE, cores = 4)
##11 seqs
forw_ali <- virtualPCR(rdb, up = SYMPSBANCRF, trimprimers = FALSE, cores = 4) #44
rev_ali <- virtualPCR(rdb, up = rc(SYMPSBANCRR), trimprimers = FALSE, cores = 4) #77

keeps <- names(forw_ali) %in% names(rev_ali) #19
forw_ali <- forw_ali[keeps]
keeps <- names(rev_ali) %in% names(forw_ali) #19
rev_ali <- rev_ali[keeps]
rev_ali <- rev_ali[names(forw_ali)]

forw_ali <- shave(forw_ali, char2dna(SYMPSBANCRF))
rev_ali <- shave(rev_ali, char2dna(rc(SYMPSBANCRR)))

concats <- paste0(sub(".+\\|", "", names(forw_ali)), dna2char(forw_ali), dna2char(rev_ali))
discards <- duplicated(concats) #10
forw_ali <- forw_ali[!discards]
rev_ali <- rev_ali[!discards]

forw_ali <- aphid::align(forw_ali)
rev_ali <- aphid::align(rev_ali)

taxids <- as.integer(sub(".+\\|", "", rownames(forw_ali)))
get_lineage(taxids, tdb)
sub("\\|.+", "", rownames(forw_ali))

tmp1 <- ape::as.character.DNAbin(forw_ali)
tmp2 <- ape::as.character.DNAbin(rev_ali)




## virtual PCR to catch sequences that include Pochon14 primer bind region
z <- virtualPCR(rdb, up = psbA_1.0, down = psbA_3.0, trimprimers = TRUE, cores = 4)

# 26 sequences retained
## derive profile hidden Markov model
model <- aphid::derivePHMM(z, cores = 4) ## 670 modules
## virtualFISH to include those without primer bind regions (i.e. including Pochon's sequence)
saveRDS(model, file = "model.rds")

## check for high scoring sequences among those initially rejected
rejects <- !(names(rdb) %in% names(z))
rejects <- rdb[rejects]
seconds <- virtualFISH(rejects, probe = model, minscore = 50,
                       minamplen = 500L, maxamplen = 1000L, cores = 4L)
## 35 retained
z <- c(z, seconds)
taxids <- as.integer(sub(".+\\|", "", names(z)))

lineages <- get_lineage(taxids, tdb)
vapply(lineages, tail, "", 1, USE.NAMES = FALSE)

## our primers are outside range
# virtualPCR(z, up = SYMPSBANCRF, trimprimers = TRUE)
# virtualPCR(z, up = rc(SYMPSBANCRR), trimprimers = TRUE)


virtualPCR(symtrim, up = SYMPSBANCRF, trimprimers = TRUE)
virtualPCR(symtrim, up = rc(SYMPSBANCRR), trimprimers = TRUE)
aphid::Viterbi(symtrim[1], z[1], type = "semiglobal")$path
aphid::Viterbi(symtrim[1], char2dna(SYMPSBANCRF), type = "semiglobal")$path
aphid::Viterbi(symtrim[1], char2dna(rc(SYMPSBANCRR)), type = "semiglobal")$path



HG515017_CDS <- structure(sym["HG515017|431157"][[1]][1:1029], class = "DNAbin")
symtrim <- virtualFISH(sym, probe = HG515017_CDS, cores = 4) # 7seqs
