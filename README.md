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

## 4 Mapping ENSEMBL IDs to HGNC symbols

#### Preparations: packages, functions, files

For mapping our Illumina ID's to HUGO symbol's we will use the Bioconductor object **'illuminaHumanv4.db'** (documentation for the object can be read from [here](http://bioconductor.org/packages/release/data/annotation/manuals/illuminaHumanv4.db/man/illuminaHumanv4.db.pdf)) which maps Illumina ID's to its matching HUGO symbol. This package will be downloaded via BiocManager as it is also a Bioconductor object. 
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("illuminaHumanv4.db", version = "3.8")
library("illuminaHumanv4.db")
```

Using the vector of probe ID's from the data sourcing script (**'GEO_tissues_data_download.R'**), we will 

---- 

# 5 Data statistics (quantile normalization)

----

## 7 References




