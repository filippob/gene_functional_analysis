##### installazione pacchetti
# BiocManager::install(c("clusterProfiler"))

## Rscript <createEntrezIDs> --projectfolder <project base folder path> --genesymbols <sample_path> -n <column_name> --outfolder <outf> --species <species> --level <level>

library("dplyr")
library("mygene")
library("ggplot2")
library("optparse")
library("data.table")
library("clusterProfiler")


####### parsing dei parametri 

option_list <- list(
  make_option(c("-p", "--projectfolder"), action="store", type="character",
              help="project base folder"),
  make_option(c("-g", "--genesymbols"), action="store", type="character",
              help="input file (one column of gene symbols)"),
  make_option(c("-n", "--columnname"), action="store", type="character",
              help="name of column of gene symbols in the input file"),
  make_option(c("-o", "--outfolder"), action="store", default="outfiles",
              help="output csv file prefix [default=outfile]"),
  make_option(c("-s", "--species"), action="store", default='human', type="character",
              help="species (between 'Anopheles, Arabidopsis, Bovine, Worm, Canine, Fly, Zebrafish, Ecoli_strain_K12, Ecoli_strain_Sakai, Chicken, Human, Mouse, Rhesus, Malaria, Chimp, Rat, Yeast, Pig, Xenopus') [default=human]"),
  make_option(c("-l", "--level"), action="store", default=3, type="integer",
              help="set level of GO depth [default=3]")
)

parser <- OptionParser(usage="%prog [options] ", option_list=option_list)
args <- parse_args(parser, positional_arguments = TRUE)
opt <- args$options
prjfolder = opt$projectfolder
infilename = opt$genesymbols
column = opt$columnname
species = gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", opt$species, perl = TRUE) ## make sure that the first letter is uppercase
outfolder = opt$outfolder
level = opt$level
#
# prjfolder = "/home/filippo/Documents/salvo/val_belice/GWAS"
# infilename = "post_gwas/gwas_prolificity_sheep_genes.csv"
# column = "uniprot_gn_symbol"
# species = "Sheep"
# outfile = "post_gwas/"
# level = 3
# taxonomy = 9940 ## sheep


# source("https://bioconductor.org/biocLite.R")
# if the org.db of the species of interest has not been installed via biocLite(), you need to do so, otherwise the script will fail
switch(species, 
       Anopheles={library("org.Ag.eg.db")},
       Arabidopsis={library("org.At.tair.db")},
       Bovine={library("org.Bt.eg.db");orgDB="org.Bt.eg.db";taxonomy=9913},
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
       #Rice={library("MeSH.Osa.eg.db");orgDB="MeSH.Osa.eg.db";taxonomy=4530},
       Xenopus={library("org.Xl.eg.db")},
       Sheep = {taxonomy = 9940},
       {                    # altrimenti ....
         stop("Specify one of this --specie options: Anopheles, Arabidopsis, Bovine, Worm, Canine, Fly, Zebrafish, Ecoli_strain_K12, Ecoli_strain_Sakai, Chicken, Human, Mouse, Rhesus, Malaria, Chimp, Rat, Yeast, Pig, Xenopus", call.=FALSE)
       }
)


# carico dal file i genesymbols
fname = file.path(prjfolder, infilename)
df <- fread(file = fname)

geneSymbolList <- df |> select({{column}}) |> filter(if_any(everything(), ~ .x != ""))

print("List of gene symbols")
print(geneSymbolList)


##### prima devo trasformare i gensymbols in entrezid 
## provo con mygene

ge <- as.data.frame(
  queryMany(
    geneSymbolList, 
    scopes="symbol", #official gene symbol (e.g. cdk2) 
    fields="entrezgene", #Entrez gene id will be returned
    species=taxonomy) #9913=Bos taurus #9606=human
)

geneEntrezList <- ge$entrezgene[!is.na(ge$entrezgene)]

print("List of Entrez Gene Codes")
print(geneEntrezList)

outfileEntrez <- paste(species, "_EntrezIDs", sep = "")
fname = file.path(prjfolder, outfolder, outfileEntrez)

fwrite(x = as.data.frame(geneEntrezList), file = fname)

print(paste("output written to", fname))

print("DONE!")




  
