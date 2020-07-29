library(rentrez)

RAC_blast <- read.delim(file="D:/bioinformatics/wgs_unmapped_blast/RAC_unmappedBLAST.out", header=FALSE)
colnames(RAC_blast) <- c("query_acc_ver", "subject_acc_ver", "pct_identity", "aln_length",
                         "num_mismatch", "gap_openings", "query_aln_start", "query_aln_end",
                         "subject_aln_start", "subject_aln_end", "e_value", "bit_score")

# entrez_search(db="nucleotide", term=paste(RAC_blast$subject_acc.ver[2], "[ACCN]", sep=""))$id