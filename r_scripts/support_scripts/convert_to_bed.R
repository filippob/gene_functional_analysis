
## script to create bed file for liftover: https://genome.ucsc.edu/cgi-bin/hgLiftOver

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
    inpfile = "paper/significant_snps.csv",
    outfile = "paper/significant_snps.bed",
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
print(paste("input file is:", config$inpfile))
print(paste("output file is:", config$outfile))


### reading data and converting format
writeLines(" - reading snp data")
fname = file.path(config$basefolder, config$inpfile)
snps <- fread(fname)

writeLines(" - converting to bed format")
snps$Chrom <- paste("chr",snps$Chrom,sep="")
snps <- snps |> 
  dplyr::select(Chrom, Position) |>
  mutate(end = Position+1) ## add 1 bp to get end position (this is the correct way, for SNP rs429585659 checked below
                           ## https://www.ensembl.org/Ovis_aries/Variation/Explore?db=core;r=1:16765815-16766815;v=rs429585659;vdb=variation;vf=50615591)

writeLines(" - writing out converted file")
fname = file.path(config$basefolder, config$outfile)
fwrite(x = snps, file = fname, col.names = FALSE, sep = "\t")

print("DONE!")

