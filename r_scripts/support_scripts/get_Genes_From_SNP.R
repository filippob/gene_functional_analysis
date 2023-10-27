##############################################################
## SCRIPT TO READ DATA AND FIND GENES CLOSE TO SNP (FROM GWAS)
##############################################################

## SNP input file must have the following three columns in position 1, 2 and 3 (first three columns):
## 1) chromosome
## 2) SNP name
## 3) position (bps)

# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("biomaRt")
# BiocManager::valid()
# BiocManager::install(version = "3.17", update=TRUE, lib = "/home/filippo/R/x86_64-pc-linux-gnu-library/4.3")

# getNamespaceInfo("BiocManager", "path")

# library("plyr")
library("plyr")
library("dplyr")
library("RCurl")
library("biomaRt")
library("ggplot2")
#library("tidyverse")
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
    outfile = "post_gwas/gwas_prolificity_sheep_genes.csv",
    ensembl_genome = "oaries_gene_ensembl",
    biomart_db = "ensembl",
    window = 50000,
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
print(paste("biomart database:", config$biomart_db))
print(paste("ensembl dataset:", config$ensembl_genome))
print(paste("window (bps):", config$window))


# listEnsembl()
###########################
## ENSEMBL object
###########################
writeLines(' - choose databse')
# listMarts()
ensembl = biomaRt::useMart(config$biomart_db) ## choose database to use (ensembl)

writeLines(' - choose dataset')
datasets <- biomaRt::listDatasets(ensembl, verbose = TRUE) # show all the possible databases on Ensembl
print(head(datasets))

# ensembl <- useEnsembl(biomart = "ensembl", dataset = "cfamiliaris_gene_ensembl")
ensembl = biomaRt::useEnsembl(biomart=config$biomart_db, dataset=config$ensembl_genome)
print(ensembl)

####################################################
## READ DATA AND FIND GENES CLOSE TO SNP (FROM GWAS)
####################################################
writeLines(' - read SNP data')
### READ DATA ##
## text file with snp name/id, chromosome, position (bps)
# input_file_name <- getURL("https://raw.githubusercontent.com/filippob/introduction_to_gwas/master/example_data/dogs_imputed.raw_cleft_lip_GWAS.results")
# results = read.table(text=input_file_name,sep=",",header=T,colClasses = c("character","integer","integer","numeric"))
# names(results)[1] <- "SNP"

## GWAS Survival to covid-19
fname = file.path(config$basefolder, config$inpfile)
res_gwas = fread(fname)
names(res_gwas)[1] <- "CHR"
names(res_gwas)[2] <- "SNP"
names(res_gwas)[3] <- "BP"

res_gwas$CHR = gsub("chr","",res_gwas$CHR)
# res_cov$BP = str_split_i(string = res_cov$SNP, pattern = ":", 2) ## string, pattern_to_split, element to return

results <- res_gwas |>
  dplyr::select(SNP,CHR,BP) |>
  # dplyr::rename(P = pvalue) |>
  mutate(BP = as.numeric(BP))

##############################################
## IF YOU WANT TO CALCULATE ADJUSTED P-VALUES
## AND FILTER SNP ACCORDINGLY
## UNCOMMENT AND RUN THE CODE BELOW
##############################################

##Calculate FDR and filter SNPs
# results$Padj <- p.adjust(results$P, method="fdr" )
#results$Padj<-p.adjust( results$P, method="bonferroni" )

##Calculate Bonferroni correction P-value
# Bonf <- 0.05/dim(results)[1]
# print(paste("The significant p-value after Bonferroni correction is",Bonf,sep=" "))

#Filter significant SNPs based on Bonferroni
 ##results <- results[results$P<Bonf,]
#Filter significant SNPs based on FDR
# fdr <- max(results$P[which(results$Padj<0.05)])
# print(paste("The significant p-value after FDR correction is",fdr,sep=" "))
# results <- results[results$Padj < fdr,]

##############################################
##############################################

# results$P <- NULL
results <- unique(results)
row.names(results) <- results$SNP

print("Unique SNP names:")
print(row.names(results))

## listAttributes(ensembl) # show the attributes of the database
writeLines(" - available attributes (head)")
attributes <- listAttributes(ensembl)  # show the attributes of the database
print(head(attributes))

genes = list()

# The getBM() function has three arguments that need to be introduced: filters, attributes and values. 
# Filters define a restriction on the query. For example you want to restrict the output to all genes 
# located on the human X chromosome then the filter chromosome_name can be used with value ‘X’. 
# The listFilters() function shows you all available filters in the selected dataset.

# c("chromosome_name","start","end") %in% listFilters(ensembl)$name

# The getBM() function is the main query function in biomaRt. It has four main arguments:
# - attributes: is a vector of attributes that one wants to retrieve (= the output of the query).
# - filters: is a vector of filters that one will use as input to the query.
# - values: a vector of values for the filters. 
# In case multple filters are in use, the values argument requires a list of values where each position in the list 
# corresponds to the position of the filters in the filters argument (see examples below).
# mart: is an object of class Mart, which is created by the useMart() function.

# getBM(attributes=c('ensembl_gene_id', 'entrezgene_id'), 
#       filters = c('chromosome_name'), 
#       values = c(snp$CHR,snp$BP-window,snp$BP+window), 
#       mart = ensembl)

writeLines(' - get genes (ensembl attributes) from SNP positions')

i = 1
for (snp_name in rownames(results)) {
    snp = results[SNP == snp_name, , ]
    print(paste("SNP n.",i, "with name", snp_name))
    genes[[snp_name]] = biomaRt::getBM(attributes = c('ensembl_gene_id',
                                         'entrezgene_id',
                                         'external_gene_name',
                                         'start_position',
                                         'end_position',
                                         'uniprot_gn_symbol',
                                         'uniprotsptrembl'
                                         # 'uniprotswissprot'
                                         ),  
                                       filters = c("chromosome_name","start","end"),
                                       values=list(snp$CHR,snp$BP-config$window,snp$BP+config$window),
                                       mart=ensembl)
    i = i+1
}

## GET GENES
#convert list to dataframe
gwas_genes <- plyr::ldply(genes, function(x) {
    rbind.data.frame(x)
})

writeLines(' - write out results')

outdir = gsub(basename(config$outfile), "", config$outfile) 
dir.create(file.path(config$basefolder, outdir), showWarnings = FALSE)

fname = file.path(config$basefolder, config$outfile)
fwrite(x = gwas_genes, file = fname)

print("DONE!")
