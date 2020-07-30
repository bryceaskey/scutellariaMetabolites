library(rentrez)

set_entrez_key("71ed3fbd7d429ed0e6409d0d0d0c0f3b3a09")

ALT_blast <- read.delim(file="/ufrc/lee/braskey/Data/WGS/bwa/BLAST_results/ALT_unmappedBLAST.out", header=FALSE)
colnames(ALT_blast) <- c("query_acc_ver", "subject_acc_ver", "pct_identity", "aln_length",
                         "num_mismatch", "gap_openings", "query_aln_start", "query_aln_end",
                         "subject_aln_start", "subject_aln_end", "e_value", "bit_score")
ALT_blast <- ALT_blast[!duplicated(ALT_blast[,"query_acc_ver"]),]

get_sci_name <- function(subject_acc_ver){
  nt_id <- entrez_search(db="nucleotide", term=paste(subject_acc_ver, "[ACCN]", sep=""))$id
  tax_id_ncbi <- entrez_summary(db="nucleotide", id=nt_id)$taxid
  tax_id <- entrez_search(db="taxonomy", term=paste(tax_id_ncbi, "[UID]", sep=""))$id
  sci_name <- entrez_summary(db="taxonomy", id=tax_id)$scientificname
  return(sci_name)
}

ALT_blast$subject_sci_name <- lapply(ALT_blast$subject_acc_ver, get_sci_name)
saveRDS(ALT_blast, file="/ufrc/lee/braskey/Data/WGS/bwa/BLAST_results/ALT_unmappedBLAST.rds")