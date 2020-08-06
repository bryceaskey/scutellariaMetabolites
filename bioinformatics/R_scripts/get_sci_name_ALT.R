library(rentrez)

set_entrez_key("71ed3fbd7d429ed0e6409d0d0d0c0f3b3a09")

ALT_blast <- read.delim(file="D:/bioinformatics/wgs_unmapped_blast/ALT_unmappedBLAST.out", header=FALSE)
colnames(ALT_blast) <- c("query_acc_ver", "subject_acc_ver", "pct_identity", "aln_length",
                         "num_mismatch", "gap_openings", "query_aln_start", "query_aln_end",
                         "subject_aln_start", "subject_aln_end", "e_value", "bit_score")
ALT_blast <- ALT_blast[!duplicated(ALT_blast[,"query_acc_ver"]),]

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

ALT_blast$subject_tax_id_ncbi <- NA
ALT_blast$subject_sci_name <- NA

for(i in 64835:nrow(ALT_blast)){
  while(TRUE){
    ALT_blast$subject_tax_id_ncbi[i] <- get_tax_id_ncbi(ALT_blast$subject_acc_ver[i])
    if(!is(ALT_blast$subject_tax_id_ncbi[i], 'try-error')) break
  }
  while(TRUE){
    ALT_blast$subject_sci_name[i] <- get_sci_name(ALT_blast$subject_tax_id_ncbi[i])
    if(!is(ALT_blast$subject_sci_name[i], 'try-error')) break
  }
  print(i)
}

#ALT_blast$subject_tax_id_ncbi <- lapply(ALT_blast$subject_acc_ver, get_tax_id_ncbi)
#ALT_blast$subject_sci_name <- lapply(ALT_blast$subject_tax_id_ncbi, get_sci_name)
#saveRDS(ALT_blast, file="D:/bioinformatics/wgs_unmapped_blast/ALT_unmappedBLAST.rds")