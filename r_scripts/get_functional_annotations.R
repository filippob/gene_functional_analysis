# Rscript <script_name> --projectfolder <project base folder path> --gene_names <geneEntrezList> --outfolder <outfolder> --species <species> --alias <alias> 

library("AnnotationHub")
library("data.table")
library("wordcloud")
library("magrittr")
library("optparse")
library("ggplot2")
library("biomaRt")
library("stringr")

####### parsing dei parametri 

option_list <- list(
  make_option(c("-p", "--projectfolder"), action="store", type="character",
              help="project base folder"),
  make_option(c("-e", "--gene_names"), action="store", type="character",
              help="input file (one columns of gene Entrez IDs, Ensembl IDs and gene symbols)"),
  make_option(c("-o", "--outfolder"), action="store", default="outfiles",
              help="output folder (path to, relative to project folder) [default=outfiles]"),
  make_option(c("-m", "--mirror"), action="store", default="",
              help="mirror biomart databse (e.g. asia.ensembl.org) [default='']"),
  make_option(c("-s", "--species"), action="store", default='human', type="character",
              help="species (between 'Anopheles, Arabidopsis, Bovine, Worm, Canine, Fly, Zebrafish, Ecoli_strain_K12, Ecoli_strain_Sakai, Chicken, Human, Mouse, Rhesus, Malaria, Chimp, Rat, Yeast, Pig, Xenopus') [default=human]"),
  make_option(c("-l", "--alias"), action="store", default="ensembl_gene_id", type="character",
              help="choose gene alias to use")
)


parser <- OptionParser(usage="%prog [options] ", option_list=option_list)
args <- parse_args(parser, positional_arguments = TRUE)
opt <- args$options
prjfolder = opt$projectfolder
infilename = opt$gene_names
species = gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", opt$species, perl = TRUE) ## make sure that the first letter is uppercase
outfolder = opt$outfolder
gene_alias = opt$alias
mirror_db = opt$mirror # (e.g. "asia.ensembl.org")


prjfolder = "/home/filippo/Documents/salvo/val_belice/GWAS"
infilename = "post_gwas/gwas_prolificity_sheep_genes.csv"
species="Goat"
outfolder ="post_gwas"
gene_alias = "uniprot_gn_symbol"
mirror_db = "useast.ensembl.org"

## select species
print(paste("selected species is", species))
switch(species, 
       Anopheles={library("org.Ag.eg.db")},
       Arabidopsis={library("org.At.tair.db")},
       Cow = {orgDB="org.Bt.eg.db";taxonomy=9913; dataset = "btaurus_gene_ensembl"},
       Worm={library("org.Ce.eg.db");orgDB="org.Ce.eg.db";taxonomy=6239},
       Canine={library("org.Cf.eg.db")},
       Fly={library("org.Dm.eg.db")},
       Zebrafish={library("org.Dr.eg.db")},
       Ecoli_strain_K12={library("org.EcK12.eg.db")},
       Ecoli_strain_Sakai={library("org.EcSakai.eg.db")},
       Chicken={library("org.Gg.eg.db");orgDB="org.Gg.eg.db";taxonomy=9031},
       Human={library("org.Hs.eg.db");orgDB="org.Hs.eg.db";taxonomy=9606},
       Mouse={library("org.Mm.eg.db");orgDB="org.Mm.eg.db";taxonomy=10090},
       Rhesus={library("org.Mmu.eg.db")},
       Malaria={library("org.Pf.plasmo.db")},
       Chimp={library("org.Pt.eg.db")},
       Rat={library("org.Rn.eg.db");orgDB="org.Rn.eg.db";taxonomy=10116},
       Yeast={library("org.Sc.sgd.db")},
       Pig={library("org.Ss.eg.db");orgDB="org.Ss.eg.db";taxonomy=9823},
       Xenopus={library("org.Xl.eg.db")},
       Goat = {taxonomy = 9925 ;sci_name = "Capra hircus"; dataset = "chircus_gene_ensembl"},
       Sheep = {taxonomy = 9940; sci_name = "Ovis aries"; dataset = "oaries_gene_ensembl"},
       {                    # altrimenti ....
         stop("Specify one of this --specie options: Anopheles, Arabidopsis, Bovine, Worm, Canine, Fly, Zebrafish, Ecoli_strain_K12, Ecoli_strain_Sakai, Chicken, Human, Mouse, Rhesus, Malaria, Chimp, Rat, Yeast, Pig, Xenopus", call.=FALSE)
       }
)


writeLines(" - reading the data")
fname = file.path(prjfolder, infilename)
genes <- fread(file=fname)

print(paste("selected gene alias is", gene_alias))
gene_list <- genes |>
  dplyr::pull(!!gene_alias) |>
  unique()

vec <- (gene_list != "")
gene_list <- gene_list[vec]

writeLines(" - selecting mart object")
## if the ensembl site is not responsive, you can try using a mirror (see below)
# mart <- useMart(biomart = "ensembl", dataset = dataset)
# mart <- useMart(biomart = "ensembl", dataset = dataset, host="useast.ensembl.org")

# ensembl <- useEnsembl(biomart = "ensembl", dataset = "cfamiliaris_gene_ensembl")
if (str_trim(mirror_db) == "") {
  print("using default biomart ensembl database")
  mart <- useMart(biomart = "ensembl", dataset = dataset)
} else {
  print(paste("using the following mirror database:", mirror_db))
  mart = useMart(biomart = "ensembl", dataset = dataset, host=mirror_db)
}


## BIOMART
writeLines(" - querying biomaRt ...")
annotated <- getBM(
  attributes = c('ensembl_gene_id', 'external_gene_name', 'uniprot_gn_symbol', 'description'), 
  filters = gene_alias,
  values = gene_list, 
  mart = mart)

writeLines(" - writing out results")

outfile <- file.path(prjfolder, outfolder, paste(species,"_functional_annotations",sep = ""))
print(paste("writing out to file", outfile))
fwrite(x=annotated, file = outfile)

print("DONE!!")