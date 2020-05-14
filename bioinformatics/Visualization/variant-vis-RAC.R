library(vcfR)
library(ape)

species <- "RAC"
variant_path <- "/ufrc/lee/braskey/Data/WGS/bwa/"
ref_path <- "/ufrc/lee/braskey/Data/ASM577160v1/ncbi_dataset/data/GCA_005771605.1/"
output_path <- "/ufrc/lee/braskey/Figures/"

pdf(paste(output_path, species, "_variantPlots.pdf", sep=""))
for (chrNum in 1:9) {
  vcf <- read.vcfR(paste(variant_path, species, "_aln_flt_chr", chrNum, ".vcf.gz", sep=""))
  ref <- ape::read.dna(paste(ref_path, "chr", chrNum, ".fa", sep=""), format="fasta")
  chrom <- create.chromR(name=paste(species, "chromosome", chrNum), vcf=vcf, seq=ref)
  chrom <- proc.chromR(chrom, verbose=TRUE)
  chromoqc(chrom, dp.alpha=20)
}
dev.off()