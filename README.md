# BCB420.2019.GEO.tissues
#### (GEO (tissue) data annotatation of human genes)

----

## 1 About this package:

----

## 2 GEO Data

#### 2.1 Data semantics

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

It will return a vector of information stored in Bioconductor's Expression Sets class (see [Bioconductor's Expression Set Class documentation](https://www.bioconductor.org/packages/3.7/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf) for more information). 

This script also returns a character vector of the Illumina probe ID's which will be used later to map the ID's to its matching HUGO symbol.

---- 

## 4 Mapping Probe IDs to HGNC symbols

#### Preparations: packages, functions, files

For mapping our Illumina ID's to HUGO symbol's we will use the Bioconductor object **'illuminaHumanv4.db'** (documentation for the object can be read from [here](http://bioconductor.org/packages/release/data/annotation/manuals/illuminaHumanv4.db/man/illuminaHumanv4.db.pdf)) which maps Illumina ID's to its matching HUGO symbol. This package will be downloaded via BiocManager as it is also a Bioconductor object. 
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("illuminaHumanv4.db", version = "3.8")
library("illuminaHumanv4.db")
```

Using the vector of probe ID's from the data sourcing script (**'GEO_tissues_data_download.R'**), we will remove any duplicates in probeID ...

```R
probeID <- MyexprNames 
length(probeID)  #699328 probe ID's extracted from the ~20 datasets
probeID <- unique(MyexprNames)
length(probeID)  #79261 unique probe ID's
```

... and map the probeID's to HGNC symbols using the **'illuminaHumanv4.db`** ...

```R
HUGOgeneAnnot <- data.frame(mapIds(illuminaHumanv4.db, probeID, "SYMBOL","PROBEID"))
colnames(HUGOgeneAnnot) <- "hgnc_symbol"

> head(HUGOgeneAnnot)
#             hgnc_symbol
#ILMN_1343291      EEF1A1
#ILMN_1343295       GAPDH
#ILMN_1651199        <NA>
#ILMN_1651209    SLC35E2A
#ILMN_1651210        <NA>
#ILMN_1651221      EFCAB1

nrow(HUGOgeneAnnot)  #79,261

```

... which returns a database with row names as the illumina probeID's and each cell in the row corresponds to the HGNC symbol it corresponds to. 


```R
sum(is.na(HUGOgeneAnnot)) #43539 probe ID's not mapped to HGNC symbol (from total 79261) = 
```

Only ~55% of the probe Id's from the 18 datasets successfully mapped to a HGNC symbol... this is not that great so we will try to map as many of the probe ID's as possible using additional mapping techniques. 

First I will collect all the probe ID's that don't have a corresponding symbol. 

```R

```

1. Mapping to synonyms
```R

```

2. Mapping to previous
```R

```

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


