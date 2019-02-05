#R Script for BCB420 GEO (tissues) data download and cleanup
#Author: Cathy Cha 
#Database: GEO (tissues) 
#
#
#
# 1. Preparations ##############################################
################################################################

#Need bioconductor package to analyse GEO datasets
if (! require(Biobase, quietly=TRUE)) {
  if (! exists("biocLite")) {
    source("https://bioconductor.org/biocLite.R")
  }
  biocLite("Biobase")
  library(Biobase)
}

if (! require(GEOquery, quietly=TRUE)) {
  if (! exists("biocLite")) {
    source("https://bioconductor.org/biocLite.R")
  }
  biocLite("GEOquery")
  library(GEOquery)
}

# 2. Load the datasets from GEO ################################
################################################################

#Interstitial fluid flow effect on noninvasive and invasive ERBB2-positive breast cancer cells
GSE64670 <- getGEO("GSE64670", GSEMatrix =TRUE, getGPL=FALSE)
GSE64670 <- GSE64670[[1]]
GSE64670Genes <- featureNames(GSE64670)

#MicroRNA-135b overexpression effect on prostate cancer cell line: time course
GSE57820 <- getGEO("GSE57820", GSEMatrix =TRUE, getGPL=FALSE)
GSE57820 <- GSE57820[[1]]
GSE57820Genes <- featureNames(GSE57820)

#Myelodysplastic syndrome disorder: bone marrow mesenchymal stromal cells
GSE61853 <- getGEO("GSE61853", GSEMatrix =TRUE, getGPL=FALSE)
GSE61853 <- GSE61853[[1]]
GSE61853Genes <- featureNames(GSE61853)

#Estradiol effect on MCF7 breast cancer cells expressing progesterone receptor-B
GSE45643 <- getGEO("GSE45643", GSEMatrix =TRUE, getGPL=FALSE)
GSE45643 <- GSE45643[[1]]
GSE45643Genes <- featureNames(GSE45643)

#miR-205 silencing effect on prostate epithelial cell line
GSE29782 <- getGEO("GSE29782", GSEMatrix =TRUE, getGPL=FALSE)
GSE29782 <- GSE29782[[1]]
GSE29782Genes <- featureNames(GSE29782)

#Short-term fasting effect on skeletal muscle: time course
GSE55924 <- getGEO("GSE55924", GSEMatrix =TRUE, getGPL=FALSE)
GSE55924 <- GSE55924[[1]]
GSE55924Genes <- featureNames(GSE55924)

#Mutant HLA-F-adjacent transcript 10 overexpression effect on HCT116 colon cancer cell line
GSE54167 <- getGEO("GSE54167", GSEMatrix =TRUE, getGPL=FALSE)
GSE54167 <- GSE54167[[1]]
GSE54167Genes <- featureNames(GSE54167)

#Helicobacter pylori infection: corpus stomach biopsies
GSE27411 <- getGEO("GSE27411", GSEMatrix =TRUE, getGPL=FALSE)
GSE27411 <- GSE27411[[1]]
GSE27411Genes <- featureNames(GSE27411)

#Rho kinase inhibition effect on epidermal keratinocyte in vitro
GSE52515 <- getGEO("GSE52515", GSEMatrix =TRUE, getGPL=FALSE)
GSE52515 <- GSE52515[[1]]
GSE52515Genes <- featureNames(GSE52515)

#Matrigel 3D culture model of JIMT1 breast cancer cells
GSE42529 <- getGEO("GSE42529", GSEMatrix =TRUE, getGPL=FALSE)
GSE42529 <- GSE42529[[1]]
GSE42529Genes <- featureNames(GSE42529)

#	miR-205 expression effect on prostate cancer cell line
GSE11701 <- getGEO("GSE11701", GSEMatrix =TRUE, getGPL=FALSE)
GSE11701 <- GSE11701[[1]]
GSE11701Genes <- featureNames(GSE11701)

#DNA methyltransferase inhibitor 5-aza-2’-deoxycytidine effect on oral cancer cell lines
GSE38823 <- getGEO("GSE38823", GSEMatrix =TRUE, getGPL=FALSE)
GSE38823 <- GSE38823[[1]]
GSE38823Genes <- featureNames(GSE38823)

#Akt inhibitor MK2206 effect on influenza H1N1 infection of non-small cell lung cancer line
GSE54293 <- getGEO("GSE54293", GSEMatrix =TRUE, getGPL=FALSE)
GSE54293 <- GSE54293[[1]]
GSE54293Genes <- featureNames(GSE54293)

#miR-542-3p overexpression effect on TP53 wild-type osteosarcoma cell line
GSE47363 <- getGEO("GSE47363", GSEMatrix =TRUE, getGPL=FALSE)
GSE47363 <- GSE47363[[1]]
GSE47363Genes <- featureNames(GSE47363)

#Dendritic cell subsets
GSE60317 <- getGEO("GSE60317", GSEMatrix =TRUE, getGPL=FALSE)
GSE60317 <- GSE60317[[1]]
GSE60317Genes <- featureNames(GSE60317)

#Long noncoding MALAT1 RNA deficiency effect on normal diploid fibroblasts W138
GSE44240 <- getGEO("GSE44240", GSEMatrix =TRUE, getGPL=FALSE)
GSE44240 <- GSE44240[[1]]
GSE44240Genes <- featureNames(GSE44240)

#Ankylosing spondylitis: blood
GSE25101 <- getGEO("GSE25101", GSEMatrix =TRUE, getGPL=FALSE)
GSE25101 <- GSE25101[[1]]
GSE25101Genes <- featureNames(GSE25101)

#3D-culture effect on fallopian tube secretory epithelial cells
GSE51220 <- getGEO("GSE51220", GSEMatrix =TRUE, getGPL=FALSE)
GSE51220 <- GSE51220[[1]]
GSE51220Genes <- featureNames(GSE51220)

#put all probe ID's into a vector 

MyexprNames <- c(GSE64670Genes, GSE57820Genes, GSE61853Genes, GSE45643Genes, GSE29782Genes, 
                 GSE55924Genes, GSE54167Genes, GSE27411Genes, GSE52515Genes, GSE42529Genes, 
                 GSE11701Genes, GSE38823Genes, GSE54293Genes, GSE47363Genes, GSE60317Genes, 
                 GSE44240Genes,GSE25101Genes, GSE51220Genes)