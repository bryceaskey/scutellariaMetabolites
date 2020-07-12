library(rentrez)

RAC_blast <- read.delim(file="D:/bioinformatics/RAC_unmappedBLAST.out", header=FALSE)
colnames(RAC_blast) <- c("query_acc.ver", "subject_acc.ver", "pct_identity", "aln_length",
                         "num_mismatch", "gap_openings", "query_aln_start", "query_aln_end",
                         "subject_aln_start", "subject_aln_end", "e_value", "bit_score")