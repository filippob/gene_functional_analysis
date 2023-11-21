##### installazione pacchetti
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("clusterProfiler")

# Rscript <genesymbols_to_GO_clusterprofiler> --projectfolder <project base folder path> --entrez_ids <geneEntrezList> --outfolder <outfolder> --species <species> --level <level> --color <variable_to_color_plots>
  
library("clusterProfiler") ## v.3.18.0
library("AnnotationHub")
library("data.table")
library("wordcloud")
library("doseplot")
library("magrittr")
library("optparse")
library("ggplot2")
library("mygene")
library("topGO")

####### parsing dei parametri 

option_list <- list(
  make_option(c("-p", "--projectfolder"), action="store", type="character",
              help="project base folder"),
  make_option(c("-e", "--entrez_ids"), action="store", type="character",
              help="input file (one column of gene Entrez IDs)"),
  make_option(c("-c", "--color"), action="store", type="character", default="pvalue",
              help="variable to use for plots (pvalue [default], p.adjust, qvalue)"),
  make_option(c("-o", "--outfolder"), action="store", default="outfiles",
              help="output folder (path to, relative to project folder) [default=outfiles]"),
  make_option(c("-s", "--species"), action="store", default='human', type="character",
              help="species (between 'Anopheles, Arabidopsis, Bovine, Worm, Canine, Fly, Zebrafish, Ecoli_strain_K12, Ecoli_strain_Sakai, Chicken, Human, Mouse, Rhesus, Malaria, Chimp, Rat, Yeast, Pig, Xenopus') [default=human]"),
  make_option(c("-l", "--level"), action="store", default=3, type="integer",
              help="set level of GO depth [default=3]")
)


parser <- OptionParser(usage="%prog [options] ", option_list=option_list)
args <- parse_args(parser, positional_arguments = TRUE)
opt <- args$options
prjfolder = opt$projectfolder
infilename = opt$entrez_ids
species = gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", opt$species, perl = TRUE) ## make sure that the first letter is uppercase
outfolder = opt$outfolder
level = opt$level
color = opt$color

# prjfolder = "/home/filippo/Documents/salvo/val_belice/GWAS"
# infilename = "post_gwas/Sheep_EntrezIDs"
# species="Sheep"
# outfolder ="post_gwas"
# level=3
# color="pvalue" # "p.adjust", "pvalue", "qvalue"


# source("https://bioconductor.org/biocLite.R")
# if the org.db of the species of interest has not been installed via biocLite(), you need to do so, otherwise the script will fail
print(paste("selected species is", species))
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
       Xenopus={library("org.Xl.eg.db")},
       enopus={library("org.Xl.eg.db")},
       Goat = {taxonomy = 9925 ;sci_name = "Capra hircus"; dataset = "chircus_gene_ensembl"},
       Sheep = {taxonomy = 9940; sci_name = "Ovis aries"},
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

### create org.db if not available
available_species = c("Anopheles", "Arabidopsis", "Bovine", "Worm", "Canine", "Fly", "Zebrafish", "Ecoli_strain_K12", "Ecoli_strain_Sakai", 
"Chicken", "Human", "Mouse", "Rhesus", "Malaria", "Chimp", "Rat", "Yeast", "Pig", "Xenopus")

print("selecting org.db")
if (!(species %in% available_species)) {
  
  hub <- AnnotationHub()
  temp = AnnotationHub::query(hub, sci_name)
  
  sci_name = gsub(" ","_",sci_name)
  db_name = paste("org.", sci_name,".eg.sqlite", sep="") 
  
  vec <- temp$title %in% db_name
  tmp <- temp[vec]
  
  org_db = hub[[names(tmp)]]
  # Oaries <- ah[["AH111977"]]
} else { org_db = orgDB}

print(org_db)

############## ora faccio la GO classification ###########

#GO Enrichment Analysis of a gene set. Given a vector of genes, 
#this function will return the enrichment GO categories after FDR control. <--ONLY ENTREZID-->
#CC=cellular component || BP= biological process || MF=molecular function 
#pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"

print("GO terms - CC")
# calcolo ed esporto i Cellular Component
ggoCC <- groupGO(gene     = geneEntrezList,
                 OrgDb    = org_db, # calls the annotation db
                 ont      = "CC",      # cellular component
                 level    = 3,         # level requested
                 readable = TRUE)

## The Gene Ratio indicates the number of genes in the input list which refer to each specific GO term
ggoCC_sort <- ggoCC@result[order(ggoCC@result$Count, decreasing=TRUE),]
outfileCC <- file.path(prjfolder, outfolder, paste(species,"_GO.CC",sep = ""))
fwrite(x = ggoCC_sort, file = outfileCC)

# barplot of GO terms
fileName <- file.path(prjfolder, outfolder,paste(species, "_GO.CC_barplot.png", sep=""))
png(fileName)
barplot(ggoCC, drop=TRUE, showCategory = 12, order = TRUE)
dev.off()

# wordcloud of GO terms
fileName <- file.path(prjfolder, outfolder, paste(species, "_GO.CC_wordcloud.png", sep=""))
png(fileName)
wordcloud(words = ggoCC@result$Description, freq = ggoCC@result$Count, min.freq = 1)
dev.off()

# calcolo ed esporto i Biological Processes
ggoBP <- groupGO(gene     = geneEntrezList,
                 OrgDb    = org_db, # calls the annotation db
                 ont      = "BP",      # cellular component
                 level    = level,         # level requested
                 readable = TRUE)

ggoBP_sort <- ggoBP@result[order(ggoBP@result$Count,decreasing=TRUE),]
outfileBP <- file.path(prjfolder, outfolder, paste(species,"_GO.BP",sep = ""))
fwrite(x=ggoBP_sort, file = outfileBP)

fileName <- file.path(prjfolder, outfolder, paste(species, "_GO.BP_barplot.png", sep=""))
png(fileName)
barplot(ggoBP, drop=TRUE, showCategory=12, order = TRUE)
dev.off()

# wordcloud of GO terms
fileName <- file.path(prjfolder, outfolder, paste(species, "_GO.BP_wordcloud.png", sep=""))
png(fileName)
wordcloud(words = ggoBP@result$Description, freq = ggoBP@result$Count, min.freq = 1)
dev.off()

# calcolo ed esporto i Molecular Function
ggoMF <- groupGO(gene     = geneEntrezList,
                 OrgDb    = org_db, # calls the annotation db
                 ont      = "MF",      # cellular component
                 level    = level,         # level requested
                 readable = TRUE)

ggoMF_sort <- ggoMF@result[order(ggoMF@result$Count,decreasing=TRUE),]
outfileMF <- file.path(prjfolder, outfolder, paste(species,"_GO.MF",sep = ""))
fwrite(x=ggoMF_sort, file = outfileMF)

fileName <- file.path(prjfolder, outfolder, paste(species, "_GO.MF_barplot.png", sep=""))
png(fileName)
barplot(ggoMF,drop=TRUE, showCategory=12, order = TRUE)
dev.off()

# wordcloud of GO terms
fileName <- file.path(prjfolder, outfolder, paste(species, "_GO.MF_wordcloud.png", sep=""))
png(fileName)
wordcloud(words = ggoMF@result$Description, freq = ggoMF@result$Count, min.freq = 1)
dev.off()

############### ora faccio la GO enrichment ##########

# calcolo ed esporto i Cellular Component
egoCC <- enrichGO(
  gene = geneEntrezList,
  OrgDb = org_db,      # calls the annotation db
  pvalueCutoff = 0.50,     
  qvalueCutoff = 0.50, 
  ont = "CC",           
  readable = TRUE
)

egoCC_sort <- egoCC@result[order(egoCC@result$Count,decreasing=TRUE),]
outfile_egoCC <- file.path(prjfolder, outfolder, paste(species,"_eGO.CC",sep = ""))
fwrite(x=egoCC_sort, file = outfile_egoCC)

fileName <- file.path(prjfolder, outfolder, paste(species, "_eGO.CC_barplot.png", sep=""))
png(fileName, height = 800)
barplot(egoCC, drop=TRUE, order=TRUE, col=color, showCategory = 10)
dev.off()

fileName <- file.path(prjfolder, outfolder, paste(species, "_eGO.CC_dotplot.png", sep=""))
png(fileName)
dotplot(object = egoCC, x = "GeneRatio", color=color)
dev.off()

# calcolo ed esporto i Biological Processes
egoBP <- enrichGO(
  gene=geneEntrezList,
  OrgDb=org_db,      # calls the annotation db
  pvalueCutoff=0.5,     
  qvalueCutoff=0.75, 
  ont = "BP",           
  readable=TRUE
)

egoBP_sort <- egoBP@result[order(egoBP@result$Count,decreasing=TRUE),]
outfile_egoBP <- file.path(prjfolder, outfolder, paste(species,"_eGO.BP",sep = ""))
fwrite(x=egoBP_sort, file = outfile_egoBP)

fileName <- file.path(prjfolder, outfolder, paste(species, "_eGO.BP_barplot.png", sep=""))
png(fileName)
barplot(egoBP, drop=TRUE, showCategory=10, col=color)
dev.off()

fileName <- file.path(prjfolder, outfolder, paste(species, "_eGO.BP_dotplot.png", sep=""))
png(fileName)
dotplot(egoBP,color=color)
dev.off()

# calcolo ed esporto i Molecular Function
egoMF <- enrichGO(
  gene=geneEntrezList,
  OrgDb=org_db,      # calls the annotation db
  pvalueCutoff=0.5,     
  qvalueCutoff=0.75, 
  ont = "MF",           
  readable=TRUE
)

egoMF_sort <- egoMF@result[order(egoMF@result$Count,decreasing=TRUE),]
outfile_egoMF <- file.path(prjfolder, outfolder, paste(species,"_eGO.MF",sep = ""))
fwrite(x=egoMF_sort, file = outfile_egoMF)

fileName <- file.path(prjfolder, outfolder, paste(species, "_GO.MF_barplot.png", sep=""))
png(fileName)
barplot(egoMF,drop=TRUE,showCategory=12, col=color)
dev.off()

fileName <- file.path(prjfolder, outfolder, paste(species, "_eGO.MF_dotplot.png", sep=""))
png(fileName)
dotplot(egoMF,color=color)
dev.off()

fileName <- file.path(prjfolder, outfolder, paste(species, "_eGO.MF_map.png", sep=""))
png(fileName)
doseplot::enrichMap(egoMF,n=25,vertex.label.font = 0.25)
dev.off()

# emapplot(egoMF)

fileName <- file.path(prjfolder, outfolder, paste(species, "_eGO.MF_graph.png", sep=""))
png(fileName)
plotGOgraph(egoMF)
dev.off()

fileName <- file.path(prjfolder, outfolder, paste(species, "_eGO.CC_graph.png", sep=""))
png(fileName)
plotGOgraph(egoCC)
dev.off()

fileName <- file.path(prjfolder, outfolder, paste(species, "_eGO.BP_graph.png", sep=""))
png(fileName)
plotGOgraph(egoBP)
dev.off()

############### END ###############


## todo: aggiungere la possibilità di stampare i plot

print("DONE!")


  
