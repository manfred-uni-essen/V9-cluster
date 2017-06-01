# script for V9 aggregation of abundance dataframes; best: following swarm_reform.R
# author: Manfred Jensen, Biodiversity group, Universität Duisburg-Essen, Germany
# V2.3  June 1, 2017   email: manfred.jensen@uni-due.de
# tested length of V9 = 120-165 BP, requested length can be modified in variable v9_end
# op_sys windows or linux, working directory can be modified; here: data in M:/

# chunk1
library(compiler)
library(data.table)
library(dplyr)
library(seqinr)
# end chunk1

# chunk2 variables
op_sys <- "windows"
filenamein <- "JAU-1_Tablesum_swo_small"
filenameout <- paste(filenamein,"_agg.csv",sep="")
filenamesorted <- paste(filenamein,"_sort.csv",sep="")
filenamefasta <- paste(filenamein,"_sort.fasta",sep="")
v9_end        <- 150
German_Excel  <- T
take_longSeqs <- F    # TRUE for longest instead of most frequent sequences to survive
wdir <- "jensen/rstudio/europa"
# end chunk2 varaibels


# chunk3 data.entry etc , linux: path to be modified
data.entry(op_sys,wdir,German_Excel)
data.entry(filenamein,v9_end,take_longSeqs)
filenameout    <- paste(filenamein,"_agg.csv",sep="")
filenamesorted <- paste(filenamein,"_sort.csv",sep="")
filenamefasta  <- paste(filenamein,"_sort.fasta",sep="")
filename       <- paste(filenamein,".csv",sep="")
T1 <- Sys.time()
if (op_sys=="windows") dbdir_M <- "M:/" else
  dbdir_M <- "/run/user/1001/gvfs/smb-share:server=132.252.96.251,share=biodiversität/"
wdir <- paste(dbdir_M,wdir,sep="")
setwd(wdir)
# end chunk3


# chunk4 read data
spe <- fread(filename,data.table=F)
colnames(spe)[1] <- "seqs"             # sequences in column 1, site abunds in following cols
# end chunk4 read data


# chunk5: sort file according to equal SSU-V9, keeping most abundant OTUs at the beginning
# ****************************************************************************************
seqs_v9end <- substring(spe$seqs,first=1,last=v9_end)
spe_ind <- as.data.table(cbind(seqs_v9end,seq(1,nrow(spe))))
names(spe_ind) <- c("seqs_v9end","rownumbers")
spe_indlist <- unstack(spe_ind,rownumbers~seqs_v9end)
spe_indlist <- unname(spe_indlist)
spe_indlist <- sapply(spe_indlist,function(x) x <- as.integer(x))
spe_inds <- sapply(spe_indlist,function(x) x[1])
spe_indlist <- spe_indlist[order(spe_inds,decreasing=F)]
inds <- unlist(spe_indlist)
spe_sorted <- spe[inds,]
v9 <- lapply(spe_indlist, function(x) sapply(x,function(y) y<-x[1]))
v9group <- unlist(v9)
seqs <- spe_sorted$seqs
spe_sorted <- select(spe_sorted,-seqs)
readsums <- rowSums(spe_sorted)
seqid <- paste("N",rownames(spe_sorted),sep="")
spe_sorted <- cbind(seqid,v9group,readsums,seqs,spe_sorted)
cat("writing  SSU_V9-sorted file\n")
filenamesorted <- paste(filenamein,"_v9sort",v9_end,".csv",sep="")

if (German_Excel) fwrite(spe_sorted, file=filenamesorted,sep=";",dec=",") else
                  fwrite(spe_sorted, file=filenamesorted,sep="\t",dec=".")
# end chunk5: preparing and writing sorted file
# ***************************************************************************************



# chunk 6: prepare fasta file from v9 sorted file
# ******************************************************************************************
seq1  <- spe_sorted$seqs  # retain seq values, prepare fasta file with seq1
seq1  <- as.list(t(seq1))    # structure needed to write fasta file
seqid <- paste(spe_sorted$seqid,"_",v9group,sep="")         # integrate group info in seqid
filenamefasta <- paste(filenamein,"_v9sort",v9_end,".fasta",sep="")    # new name
cat("writing  SSU_V9-sorted fasta file\n")                    # info
write.fasta(sequences=seq1,names=seqid,file.out=filenamefasta,open="w")
# end chunk6: writing fasta file************************************************************




# chunk7: aggregate OTUs with identical V9
# ******************************************************************************************
cat("assemble SSU_V9-aggregated file\n")
spe_v9 <- spe_sorted$v9group                      # v9 group information vector
spe_strs <- select(spe_sorted,seqid:seqs)         # extract non integer, blast infos
spenew <- select(spe_sorted,-c(seqid:seqs))       # pure abundance matrix
# aggregation
spenew$n_variants <- rep(1,nrow(spe))             # ones, to be summed up by aggregate
spenew  <-  aggregate(spenew,by=list(spe_v9),sum) # incl. n_variants, i.e. sum(n*1)
# main aggregation finished
n_variants <- spenew$n_variants                   # n_variants automatically by aggregate
spenew <- select(spenew,-n_variants,-Group.1)     # pure abundance integer matrix 
seqlength <- nchar(as.character(spe_strs$seqs))
if (take_longSeqs) 
  freq_or_long <- seqlength else 
  freq_or_long <- spe_strs$readsums
ind_min_evalue <- tapply(freq_or_long, spe_strs$v9group, which.max) # subindices within v9 group list
spe_strs_inds <- cumsum(c(0,n_variants[1:(length(n_variants)-1)])) + ind_min_evalue #determine indices
spe_strs_subs <- spe_strs[spe_strs_inds,]                    # select representative for each v9 group
seqlength <- seqlength[c(spe_strs_inds)]
readsums <- rowSums(spenew)                                      # how large is each v9 group ?
seqid <- spe_strs_subs$seqid
spe_strs_subs <- select(spe_strs_subs,-seqid)                    # seqid already stored seqid
v9group <- spe_strs_subs$v9group
seqs <- spe_strs_subs$seqs
spe_strs_subs <- cbind(seqid,v9group,n_variants,readsums,seqlength,seqs,spenew) 
# end chunk7                         # bind info and abundance matrices


# chunk8: write aggregated dataframe to file
if (take_longSeqs) filenameaggv9 <- paste(filenamein,"_v9agglongest",v9_end,".csv",sep="") else
         filenameaggv9 <- paste(filenamein,"_v9agg",v9_end,".csv",sep="")  # new name
cat("writing  SSU_V9-aggregated file: ",filenameaggv9,"\n")                # info
if (German_Excel) fwrite(spe_strs_subs,file=filenameaggv9,sep=";",dec=",") else
                  fwrite(spe_strs_subs,file=filenameaggv9,sep="\t",dec=".")
T3 <- Sys.time()
cat("Time elapsed (total ): ",T3-T1,"\n")
# end chunk8************************************************************************************


# THE END