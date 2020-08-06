library(rentrez)
library(purrr)

set_entrez_key("1f49c549b76a57b48db8b6eca1c3200afd09")

BAR_blast <- read.delim(file="D:/bioinformatics/wgs_unmapped_blast/BAR_unmappedBLAST.out", header=FALSE)
colnames(BAR_blast) <- c("query_acc_ver", "subject_acc_ver", "pct_identity", "aln_length",
                         "num_mismatch", "gap_openings", "query_aln_start", "query_aln_end",
                         "subject_aln_start", "subject_aln_end", "e_value", "bit_score")
BAR_blast <- BAR_blast[!duplicated(BAR_blast[,"query_acc_ver"]),]

get_tax_id_ncbi <- function(subject_acc_ver){
  nt_id <- entrez_search(db="nucleotide", term=paste(subject_acc_ver, "[ACCN]", sep=""))$id
  if(length(nt_id)==0){
    return(NA)
  }else{
    tax_id_ncbi <- entrez_summary(db="nucleotide", id=nt_id, always_return_list=TRUE)[[1]]$taxid
    return(tax_id_ncbi)
  }
}
get_tax_id_ncbi2 <- possibly(get_tax_id_ncbi, otherwise=NA)

get_sci_name <- function(tax_id_ncbi){
  if(is.na(tax_id_ncbi)){
    return(NA)
  }else{
    tax_id <- entrez_search(db="taxonomy", term=paste(tax_id_ncbi, "[UID]", sep=""))$id
    sci_name <- entrez_summary(db="taxonomy", id=tax_id, always_return_list=TRUE)[[1]]$scientificname
    return(sci_name)
  }
}
get_sci_name2 <- possibly(get_sci_name, otherwise=NA)

#BAR_blast$subject_tax_id_ncbi <- NA
#BAR_blast$subject_sci_name <- NA

for(i in 66673:nrow(BAR_blast)){
  print(i)
  BAR_blast$subject_tax_id_ncbi[i] <- get_tax_id_ncbi2(BAR_blast$subject_acc_ver[i])
  BAR_blast$subject_sci_name[i] <- get_sci_name2(BAR_blast$subject_tax_id_ncbi[i])
  fail_count <- 0
  while(is.na(BAR_blast$subject_tax_id_ncbi[i]) || is.na(BAR_blast$subject_sci_name[i])){
    print("Error encountered, trying again")
    Sys.sleep(5)
    BAR_blast$subject_tax_id_ncbi[i] <- get_tax_id_ncbi2(BAR_blast$subject_acc_ver[i])
    BAR_blast$subject_sci_name[i] <- get_sci_name2(BAR_blast$subject_tax_id_ncbi[i])
    fail_count <- fail_count + 1
    if(fail_count > 10){
      print(paste("After 10 tries, could not access data for:", BAR_blast$subject_acc_ver[i]))
      print("Moving to next query")
      break
    }
  }
}

fail_count <- 0
while(sum(is.na(BAR_blast$subject_sci_name)) > 0){
  if(fail_count > 0){
    Sys.sleep(500)
  }
  failed_queries <- which(is.na(BAR_blast$subject_sci_name))
  for(i in failed_queries){
    BAR_blast$subject_tax_id_ncbi[i] <- get_tax_id_ncbi2(BAR_blast$subject_acc_ver[i])
    BAR_blast$subject_sci_name[i] <- get_sci_name2(BAR_blast$subject_tax_id_ncbi[i])
  }
  fail_count <- fail_count + 1
  print(paste("Iteration", fail_count, "completed, could not access data for", sum(is.na(BAR_blast$subject_sci_name)), "accessions."))
  if(fail_count > 10){
    print(paste("After 10 tries, could not access data for:", sum(is.na(BAR_blast$subject_sci_name)), "accessions."))
    break
  }
}

#BAR_blast$subject_tax_id_ncbi <- lapply(BAR_blast$subject_acc_ver, get_tax_id_ncbi)
#BAR_blast$subject_sci_name <- lapply(BAR_blast$subject_tax_id_ncbi, get_sci_name)
#saveRDS(BAR_blast, file="/blue/lee/braskey/Data/WGS/bwa/BLAST_results/BAR_unmappedBLAST.rds")