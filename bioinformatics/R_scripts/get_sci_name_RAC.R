library(rentrez)

#set_entrez_key("afbfd9477886e9579b353069d8ca2c61cd08")

RAC_blast <- read.delim(file="/blue/lee/braskey/Data/WGS/bwa/BLAST_results/RAC_unmappedBLAST.out", header=FALSE)
colnames(RAC_blast) <- c("query_acc_ver", "subject_acc_ver", "pct_identity", "aln_length",
                         "num_mismatch", "gap_openings", "query_aln_start", "query_aln_end",
                         "subject_aln_start", "subject_aln_end", "e_value", "bit_score")
RAC_blast <- RAC_blast[!duplicated(RAC_blast[,"query_acc_ver"]),]

get_tax_id_ncbi <- function(subject_acc_ver){
  nt_id <- entrez_search(db="nucleotide", term=paste(subject_acc_ver, "[ACCN]", sep=""))$id
  if(length(nt_id)==0){
    return(NA)
  }else{
    tax_id_ncbi <- entrez_summary(db="nucleotide", id=nt_id, always_return_list=TRUE)[[1]]$taxid
    return(tax_id_ncbi)
  }
}

get_sci_name <- function(tax_id_ncbi){
  if(is.na(tax_id_ncbi)){
    return(NA)
  }else{
    tax_id <- entrez_search(db="taxonomy", term=paste(tax_id_ncbi, "[UID]", sep=""))$id
    sci_name <- entrez_summary(db="taxonomy", id=tax_id, always_return_list=TRUE)[[1]]$scientificname
    return(sci_name)
  }
}

RAC_blast$subject_tax_id_ncbi <- NA
RAC_blast$subject_sci_name <- NA

for(i in 10:nrow(RAC_blast)){
  while(TRUE){
    RAC_blast$subject_tax_id_ncbi[i] <- get_tax_id_ncbi(RAC_blast$subject_acc_ver[i])
    if(!is(RAC_blast$subject_tax_id_ncbi[i], 'try-error')) break
  }
  while(TRUE){
    RAC_blast$subject_sci_name[i] <- get_sci_name(RAC_blast$subject_tax_id_ncbi[i])
    if(!is(RAC_blast$subject_sci_name[i], 'try-error')) break
  }
  print(i)
}

RAC_blast$subject_tax_id_ncbi <- lapply(RAC_blast$subject_acc_ver, get_tax_id_ncbi)
RAC_blast$subject_sci_name <- lapply(RAC_blast$subject_tax_id_ncbi, get_sci_name)
saveRDS(RAC_blast, file="/blue/lee/braskey/Data/WGS/bwa/BLAST_results/RAC_unmappedBLAST.rds")