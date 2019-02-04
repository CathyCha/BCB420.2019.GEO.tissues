# BCB420.2019.GEO.tissues
#### (GEO (tissue) data annotatation of human genes)

----

## 1 About this package

This package will download ~20 datasets from NCBI's GEO repository and map them to HGNC symbols, quantile normalize the experiments to produce data statistics and annotate a sample gene set. 

----

## 2 GEO Data

The Gene Expression Omnibus (GEO) provided by NCBI, is a public repository for genomic data submitted by the research community. Datasets can be downloaded manually as SOFT files or MINiML files. The data download can also be downloaded in a script using the GEOquery package in Bioconductor. 

When downloaded using the getGEO function, the dataset gets saved as a Bioconductor Expression Set class (see [Bioconductor's Expression Set Class documentation](https://www.bioconductor.org/packages/3.7/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf) for more information) which includes the title of the experiment, experiment data, assay data, feature data, and many more features of the experiment conducted.

----

## 3 Data download and cleanup

#### Preparations: packages

GEOquery is a Bioconductor package that retrieves information from NCBI's Gene Expression Omnibus (GEO). It will be downloaded via biocLite since it is a bioconductor package. 

```R
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
```

The script **`Geo_tissues_data_download.R`** will include instructions on how to download the 18 tissue datasets from the GEO NCBI database using the **`getGEO()`** function from the GEOquery package.

It will return a vector of expression sets. 

This script also returns a character vector of the Illumina probe ID's which will be used later to map the ID's to its matching HUGO symbol.

---- 

## 4 Mapping Probe IDs to HGNC symbols

#### Preparations: packages, functions, files

Using the vector of probe ID's from the data sourcing script (**'GEO_tissues_data_download.R'**), we will remove any duplicates in probeID ...

```R
probeID <- MyexprNames 
length(probeID)  #699328 probe ID's extracted from the ~20 datasets
probeID <- unique(MyexprNames)
length(probeID)  #79261 unique probe ID's
```

... and map the probeID's to HGNC symbols using biomaRt...

```R
ensembl <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")

HUGOgeneAnnot <- biomaRt::getBM(attributes=c("illumina_humanht_12_v4", "hgnc_symbol", "ensembl_transcript_id"), 
      filters = "illumina_humanht_12_v4", 
      values = probeID,
      mart = ensembl)

> head(HUGOgeneAnnot)
#  illumina_humanht_12_v4 hgnc_symbol ensembl_transcript_id
#1           ILMN_1652366       NLRP7       ENST00000620820
#2           ILMN_1652366       NLRP7       ENST00000618995
#3           ILMN_1651838        RND1       ENST00000309739
#4           ILMN_1651838        RND1       ENST00000548445
#5           ILMN_1651838        RND1       ENST00000649147
#6           ILMN_1651838        RND1       ENST00000553260

nrow(HUGOgeneAnnot)  #137,830
> colnames(HUGOgeneAnnot) 
#[1] "illumina_humanht_12_v4" "hgnc_symbol"            "ensembl_transcript_id" 

#this took VERY long to map (~1hr) so make sure to save it
save(HUGOgeneAnnot, file = "HUGOMap.RData")

```

... which returns a database with columns as illumina probe id, the corresponsing HGNC symbol, and ensembl id (in order) 


```R
sum(HUGOgeneAnnot$hgnc_symbol == "")  # there are 7532 probe Id's that did not get mapped to a HGNC symbol
(sum(HUGOgeneAnnot$hgnc_symbol == "")/nrow(HUGOgeneAnnot)) *100 # 5.5% 

```

Only ~94.5% of the probe Id's from the 18 datasets successfully mapped to a HGNC symbol... to try to map the remaining 5.5% we will try alternative approaches. 


####Mapping probe ID's that don't have HGNC symbols

First I will collect all the probe ID's that don't have a corresponding symbol. 

```R
sel <- HUGOgeneAnnot$hgnc_symbol == ""
noSym <- HUGOgeneAnnot[sel,]
nrow(noSym) #7532, confirmed that we collected all empty HGNC rows
```

####1. Mapping to synonyms & previous symbols in HGNC.RData
```R
#map them to the HGNC.RData file to look for synonyms or previous symbols
sel <- noSym$ensembl_transcript_id
#see if any of the ensembl ID's missing symbols are in the HGNC dataframe
HGNCsyn <- sel %in% HGNC$EnsID
HGNCmatch <- HGNC[HGNCsyn,] #no matches
```

####2. Mapped to illuminaHumanv4.db 
An alternative method to map our Illumina ID's to HUGO symbol's using the Bioconductor object **'illuminaHumanv4.db'** (documentation for the object can be read from [here](http://bioconductor.org/packages/release/data/annotation/manuals/illuminaHumanv4.db/man/illuminaHumanv4.db.pdf)) which maps Illumina ID's to its matching HUGO symbol. This object is outdated as is used HGNC symbol data from 2015, which is why it was not used initially. 

This package will be downloaded via BiocManager as it is also a Bioconductor object. 
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("illuminaHumanv4.db", version = "3.8")
library("illuminaHumanv4.db")
```

```R
illumID <- noSym$illumina_humanht_12_v4
illumDBmatch <- data.frame(select(illuminaHumanv4.db, 
       keys = x, 
       columns=c("SYMBOL", "PROBEID"), 
       keytype="PROBEID"))

sum(is.na(illumDBmatch$SYMBOL)) #2706 (7532 - 2706 = 4826) mapped 4826 additional HGNC symbols
#1.9% not mapped = 98.1% coverage! - pretty good! 
```
Using this database we were able to map 4826 additional probe ID's to HGNC symbols. This results in a total coverage of 98.1%

The remaining 1.9% that were not mapped to an HGNC symbol could be due to (TODO: some probe id's were just not positively mapped to a gene symbol??)

####Clean up data

Now we will try to clean up the data as many of the probes id's were mapped to the same HGNC multiple times because there were multiple ensembl id's that mapped to the probe (TODO: why?)

```R


```
We are only concerned with the corresponding HGNC symbol for each probe id, so we will remove all duplicates of HGNC symbols

---- 

# 5 Data statistics (quantile normalization)

----

# 6 Annotating gene sets with GEO data 

---- 

To conclude, we annotate the example gene set, validate the annotation, and store the data

```R
#copy and pasted from STRING database and BCB420 resources 
xSet <- c("AMBRA1", "ATG14", "ATP2A1", "ATP2A2", "ATP2A3", "BECN1", "BECN2",
          "BIRC6", "BLOC1S1", "BLOC1S2", "BORCS5", "BORCS6", "BORCS7",
          "BORCS8", "CACNA1A", "CALCOCO2", "CTTN", "DCTN1", "EPG5", "GABARAP",
          "GABARAPL1", "GABARAPL2", "HDAC6", "HSPB8", "INPP5E", "IRGM",
          "KXD1", "LAMP1", "LAMP2", "LAMP3", "LAMP5", "MAP1LC3A", "MAP1LC3B",
          "MAP1LC3C", "MGRN1", "MYO1C", "MYO6", "NAPA", "NSF", "OPTN",
          "OSBPL1A", "PI4K2A", "PIK3C3", "PLEKHM1", "PSEN1", "RAB20", "RAB21",
          "RAB29", "RAB34", "RAB39A", "RAB7A", "RAB7B", "RPTOR", "RUBCN",
          "RUBCNL", "SNAP29", "SNAP47", "SNAPIN", "SPG11", "STX17", "STX6",
          "SYT7", "TARDBP", "TFEB", "TGM2", "TIFA", "TMEM175", "TOM1",
          "TPCN1", "TPCN2", "TPPP", "TXNIP", "UVRAG", "VAMP3", "VAMP7",
          "VAMP8", "VAPA", "VPS11", "VPS16", "VPS18", "VPS33A", "VPS39",
          "VPS41", "VTI1B", "YKT6")

#for our annotation we return the corresponding Illumina ID for the sample gene set
sel <- (HUGOgeneAnnot$Gene %in% xSet)
xSetProbes <- HUGOgeneAnnot[sel,]

# Statistics:
length(xSetProbes) #83 
``` 


## 7 References

Referenced the code from this Biostar's thread for part 4 (mapping to HGNC symbols)
https://www.biostars.org/p/109248/?fbclid=IwAR2kRB7vxGayHMvON4Zjvtaf_Wjw4xTwrVNte-qlwtKRzo8B-EdEhkeZkCs#109250


