##### installazione pacchetti
#source("https://bioconductor.org/biocLite.R")
#install.packages("optparse")
#biocLite("mygene")

# Rscript <genesymbols_to_Kegg_clusterprofiler> --projectfolder <project base folder path> --entrez_ids <geneEntrezList> --outfolder <outfolder> --species <species> --level <level> 

library("clusterProfiler")
library("data.table")
library("optparse")
library("ggplot2")
library("mygene")


####### parsing dei parametri 
library(optparse)
option_list <- list(
  make_option(c("-p", "--projectfolder"), action="store", type="character",
              help="project base folder"),
  make_option(c("-e", "--entrez_ids"), action="store", type="character",
              help="input file (one column of gene Entrez IDs)"),
  make_option(c("-o", "--outfolder"), action="store", default="outfolder",
              help="output csv file prefix [default=outfile]"),
  make_option(c("-s", "--species"), action="store", default='human', type="character",
              help="specie (between 'Anopheles, Arabidopsis, Bovine, Worm, Canine, Fly, Zebrafish, Ecoli_strain_K12, Ecoli_strain_Sakai, Chicken, Human, Mouse, Rhesus, Malaria, Chimp, Rat, Yeast, Pig, Xenopus') [default=human]"),
  make_option(c("-l", "--level"), action="store", default=3, type="integer",
              help="set level of GO depth [default=3]")
#  make_option(c("-n", "--count_lines"), action="store_true", default=FALSE,
#              help="Count the line numbers [default]"),
#  make_option(c("-f", "--factor"), type="integer", default=3,
#              help="Multiply output by this number [default %default]")
)

parser <- OptionParser(usage="%prog [options] ", option_list=option_list)
args <- parse_args(parser, positional_arguments = TRUE)
opt <- args$options
prjfolder = opt$projectfolder
infilename=opt$entrez_ids
species=opt$species
outfolder=opt$outfolder
level=opt$level

# prjfolder = "/home/filippo/Documents/salvo/val_belice/GWAS"
# infilename = "post_gwas/Sheep_EntrezIDs"
# species="Sheep"
# outfolder ="post_gwas"
# level=3

#source("https://bioconductor.org/biocLite.R")
print(paste("selected species is", species))
switch(species, 
       Anopheles={library("org.Ag.eg.db")},
       Arabidopsis={library("org.At.tair.db")},
       Bovine={library("org.Bt.eg.db");orgDB="org.Bt.eg.db";taxonomy=9913;keggCode="bta"},
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
       Sheep = {taxonomy = 9940; sci_name = "Ovis aries"; keggCode="oas"},
       {                    # altrimenti ....
         stop("Specify one of this --specie options: Anopheles, Arabidopsis, Bovine, Worm, Canine, Fly, Zebrafish, Ecoli_strain_K12, Ecoli_strain_Sakai, Chicken, Human, Mouse, Rhesus, Malaria, Chimp, Rat, Yeast, Pig, Xenopus", call.=FALSE)
       }
)

# carico dal file le gene Entrez IDs
fname = file.path(prjfolder, infilename)
entrezIDs <- fread(file=fname)
geneEntrezList <- unique(as.character(entrezIDs$geneEntrezList))

print("List of Entrez Gene Codes")
print(geneEntrezList)

############## ora faccio la Kegg enrichment ###########
kk <- enrichKEGG(gene         = geneEntrezList,
                 organism     = keggCode,  # calls the annotation db
                 pvalueCutoff = 0.5,
                 qvalueCutoff = 0.65)

kk_sort <- kk@result[order(kk@result$Count,decreasing=TRUE),]
outfile_kk = file.path(prjfolder, outfolder, paste(species,"_kk.csv",sep = ""))
fwrite(x=kk_sort, file = outfile_kk, sep = "\t")

############### END ###############

print("DONE!!")

## todo: aggiungere la possibilità di stampare i plot



  
