library(rentrez)
library(purrr)

set_entrez_key("1712966e23cec904480c3f55079d888c4808")

HAV_blast <- read.delim(file="D:/bioinformatics/wgs_unmapped_blast/HAV_unmappedBLAST.out", header=FALSE)
colnames(HAV_blast) <- c("query_acc_ver", "subject_acc_ver", "pct_identity", "aln_length",
                         "num_mismatch", "gap_openings", "query_aln_start", "query_aln_end",
                         "subject_aln_start", "subject_aln_end", "e_value", "bit_score")
HAV_blast <- HAV_blast[!duplicated(HAV_blast[,"query_acc_ver"]),]

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

HAV_blast$subject_tax_id_ncbi <- NA
HAV_blast$subject_sci_name <- NA

for(i in 32217:nrow(HAV_blast)){
  print(i)
  HAV_blast$subject_tax_id_ncbi[i] <- get_tax_id_ncbi2(HAV_blast$subject_acc_ver[i])
  HAV_blast$subject_sci_name[i] <- get_sci_name2(HAV_blast$subject_tax_id_ncbi[i])
  fail_count <- 0
  while(is.na(HAV_blast$subject_tax_id_ncbi[i]) || is.na(HAV_blast$subject_sci_name[i])){
    print("Error encountered, trying again")
    Sys.sleep(5)
    HAV_blast$subject_tax_id_ncbi[i] <- get_tax_id_ncbi2(HAV_blast$subject_acc_ver[i])
    HAV_blast$subject_sci_name[i] <- get_sci_name2(HAV_blast$subject_tax_id_ncbi[i])
    fail_count <- fail_count + 1
    if(fail_count > 10){
      print(paste("After 10 tries, could not access data for:", HAV_blast$subject_acc_ver[i]))
      print("Moving to next query")
      break
    }
  }
}

fail_count <- 0
while(sum(is.na(HAV_blast$subject_sci_name)) > 0){
  if(fail_count > 0){
    Sys.sleep(500)
  }
  failed_queries <- which(is.na(HAV_blast$subject_sci_name))
  for(i in failed_queries){
    HAV_blast$subject_tax_id_ncbi[i] <- get_tax_id_ncbi2(HAV_blast$subject_acc_ver[i])
    HAV_blast$subject_sci_name[i] <- get_sci_name2(HAV_blast$subject_tax_id_ncbi[i])
  }
  fail_count <- fail_count + 1
  print(paste("Iteration", fail_count, "completed, could not access data for", sum(is.na(HAV_blast$subject_sci_name)), "accessions."))
  if(fail_count > 10){
    print(paste("After 10 tries, could not access data for:", sum(is.na(HAV_blast$subject_sci_name)), "accessions."))
    break
  }
}

#HAV_blast$subject_tax_id_ncbi <- lapply(HAV_blast$subject_acc_ver, get_tax_id_ncbi)
#HAV_blast$subject_sci_name <- lapply(HAV_blast$subject_tax_id_ncbi, get_sci_name)
#saveRDS(HAV_blast, file="/blue/lee/braskey/Data/WGS/bwa/BLAST_results/HAV_unmappedBLAST.rds")