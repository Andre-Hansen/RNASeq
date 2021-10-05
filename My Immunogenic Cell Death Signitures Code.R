library(tidyverse)

#Make a list of cell death signitures in ensamble gene codes 
Immunogenic_Cell_Death_signitures <- c("ENSG00000232810", 
                                       "ENSG00000186810", 
                                       "ENSG00000089041", 
                                       "ENSG00000137752", 
                                       "ENSG00000162711", 
                                       "ENSG00000125538", 
                                       "ENSG00000154589", 
                                       "ENSG00000010610", 
                                       "ENSG00000153563", 
                                       "ENSG00000172116", 
                                       "ENSG00000180644", 
                                       "ENSG00000111537", 
                                       "ENSG00000112115", 
                                       "ENSG00000177663")

#Make a list with ID signitures 
Cell_Death_Signitures <- c("TNF", "CXCR3", "P2RX7", "CASP1", "NLRP3", "IL1B", 
                           "LY96", "CD4", "CD8A", "CD8B", "PRF1", "IFNG", "IL17A", 
                           "IL17RA")
#Testing that the matching works 
Cell_Death_Test <- c("FGR")
myListTest <- Cell_Death_Signitures
CheckTest <- subset(GeneIDTest,rownames(GeneIDTest) %in% myListTest)
dim(CheckTest)
#Make Heatmap of this 
icdmatrix<-GeneIDTest[rownames(GeneIDTest) %in% Cell_Death_Signitures,]
heatmap.2(icdmatrix)
My_Heatmap <- heatmap(icdmatrix)


#Make list 
mylist <- Immunogenic_Cell_Death_signitures
myList <- Cell_Death_Signitures
#This is to confirm the code works by using ensamble ids which are 100% in the dataset
mylist <- c("ENSG00000165092","ENSG00000270661")

#Make the list of the dataset
longlist <- rownames(vMatrix)

#Create the subset
new <- subset(vMatrix,rownames(vMatrix) %in% mylist)

#Check dimensions
dim(new)

#Another way to do the same thing
x <- which(rownames(vMatrix) %in% mylist)

new <- vMatrix[x,]

vMatrix

rownames(vMatrix) %in% Immunogenic_Cell_Death_signitures
which(rownames(vMatrix) %in% Immunogenic_Cell_Death_signitures)
vMatrix[rownames(vMatrix) %in% Immunogenic_Cell_Death_signitures,]

icdmatrix<-vMatrix[rownames(vMatrix) %in% Immunogenic_Cell_Death_signitures,]
heatmap.2(icdmatrix)
My_Heatmap <- heatmap(icdmatrix)



#Need to figure out what the heatmap means 
#Talk to mik about how to refine this down 

vMatrix = v$E

#Change the row names to samples 
n = 38
prefix <- "Sample"
suffix <- seq(1:n)
sampleNames <- paste(prefix, suffix, sep = ".")
sampleNames
colnames(vMatrix) = sampleNames

#Want to change the ensemble ID's to Gene ID's 
BiocManager::install("EnsDb.Hsapiens.v79")
library(EnsDb.Hsapiens.v79)

# 1. Convert from ensembl.gene to gene.symbol

ICDGenes <- c("ENSG00000232810", 
              "ENSG00000186810", 
              "ENSG00000089041", 
              "ENSG00000137752", 
              "ENSG00000162711", 
              "ENSG00000125538", 
              "ENSG00000154589", 
              "ENSG00000010610", 
              "ENSG00000153563", 
              "ENSG00000172116", 
              "ENSG00000180644", 
              "ENSG00000111537", 
              "ENSG00000112115", 
              "ENSG00000177663")
ICDgeneIDs <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ICDGenes, 
                              keytype = "GENEID", columns = c("SYMBOL","GENEID"))


ensembl.genes <- rownames(vMatrix)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, 
                              keytype = "GENEID", columns = c("SYMBOL","GENEID"))

# 2. Convert from gene.symbol to ensembl.gene

geneSymbols <-  c('DDX26B','CCDC83',  'MAST3', 'RPL11', 'ZDHHC20',  'LUC7L3',  
                  'SNORD49A',  'CTSH', 'ACOT8')
geneIDs2 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= geneSymbols, 
                              keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))

#Remove the ones which dont have gene ID's
GeneIDTest <- vMatrix[geneIDs1$GENEID, ]
#Then change the ensemble ID's to gene ID's
rownames(GeneIDTest) <- geneIDs1$SYMBOL


#Save Data
save.image("ICDCode.RData")
load("ICDCode.RData")
library(tidyverse)