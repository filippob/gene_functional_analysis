
library("tidyverse")
library("data.table")

args = commandArgs(trailingOnly=TRUE)
if (length(args) >= 1) {
  
  #loading the parameters
  source(args[1])
  # source("Analysis/hrr/config.R")
  
} else {
  #this is the default configuration, used for development and debug
  writeLines('Using default config')
  
  #this dataframe should be always present in config files, and declared
  #as follows
  config = NULL
  config = rbind(config, data.frame(
    #base_folder = '~/Documents/SMARTER/Analysis/hrr/',
    #genotypes = "Analysis/hrr/goat_thin.ped",
    basefolder = "/home/filippo/Documents/salvo/val_belice/GWAS",
    inpbed = "paper/hglft_genome_1e659_cefc70.bed", ## new SNP coordinates obtained from liftover: https://genome.ucsc.edu/cgi-bin/hgLiftOver
    inpold = "paper/significant_snps.csv",
    outfile = "paper/significant_snps.updated.csv",
    species = 'sheep',
    force_overwrite = FALSE
  ))
}

#############
# PARAMETERS
#############
writeLines(' - current values of parameters')
print(paste("species is:", config$species))
print(paste("base folder is:", config$basefolder))
print(paste("input bed file is:", config$inpbed))
print(paste("input snp file is:", config$inpold))
print(paste("output file is:", config$outfile))


### reading data and converting format
writeLines(" - reading snp data")
fname = file.path(config$basefolder, config$inpold)
snps <- fread(fname)

writeLines(" - reading bed file")
fname = file.path(config$basefolder, config$inpbed)
bed <- fread(fname)

## updating SNP position
writeLines(" - updating SNP position")
snps$Position <- bed$V2

writeLines(" - writing out updated file")
fname = file.path(config$basefolder, config$outfile)
fwrite(x = snps, file = fname)

print("DONE!")

