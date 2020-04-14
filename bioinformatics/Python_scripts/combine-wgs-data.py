# Script to combine WGS data into a single file
fp = open("/ufrc/lee/braskey/Data/ASM577160v1/ncbi_dataset/data/GCA_005771605.1/wgs-merged.fna", "w")

for i in range(1, 10):
  if i > 1:
    fp.write("\n")

  filename = "/ufrc/lee/braskey/Data/ASM577160v1/ncbi_dataset/data/GCA_005771605.1/chr" + str(i) + ".fna"
  
  first_line = True
  for line in open(filename, 'r'):
    if first_line == True:
      fp.write(">chr" + str(i) + "\n")
      first_line = False
    else:
      fp.write(line)

fp.close()
