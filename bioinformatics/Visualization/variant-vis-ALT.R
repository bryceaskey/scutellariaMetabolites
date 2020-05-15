library(vcfR)
library(ape)

species <- "ALT"
variant_path <- "/ufrc/lee/braskey/Data/WGS/bwa/"
ref_path <- "/ufrc/lee/braskey/Data/ASM577160v1/ncbi_dataset/data/GCA_005771605.1/"
output_path <- "/ufrc/lee/braskey/Figures/"

for (chrNum in 1:9) {
  vcf <- read.vcfR(paste(variant_path, species, "_aln_flt_chr", chrNum, ".vcf.gz", sep=""))
  ref <- ape::read.dna(paste(ref_path, "chr", chrNum, ".fa", sep=""), format="fasta")
  chrom <- create.chromR(name=paste(species, "chromosome", chrNum), vcf=vcf, seq=ref)
  chrom <- proc.chromR(chrom, verbose=TRUE)
 
  png(paste(output_path, species, "_chr", chrNum, "_variantPlot.png", sep=""),
      height=2000, width=3000)
  chromoqc(chrom, dp.alpha=150)
  dev.off()
}