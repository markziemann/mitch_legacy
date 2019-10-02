# docker run -d -p 8787:8787 -e USER=rstudio -e PASSWORD=mitch -e ROOT=TRUE rocker/rstudio
FROM rocker/rstudio

RUN apt update && \
  apt install -y \
    zlib1g-dev \
    libxml2-dev \
    tcl8.6-dev \
    tk8.6-dev \

# cran install
install.packages("devtools")
install.packages("plyr")
install.packages("tidyverse")
install.packages("reshape2")
install.packages("ashr")
install.packages("parallel")
install.packages("grid")
install.packages("gridExtra")
install.packages("beeswarm")
install.packages("GGally")
install.packages("knitr")
install.packages("markdown")
install.packages("gplots")
install.packages("purrr")
install.packages("Seurat")
install.packages("Rmisc")
install.packages("UpSetR")
install.packages("vioplot")
install.packages("rmarkdown")
install.packages("xfun")

                 

# github install
devtools::install_github("markziemann/dee2/getDEE2")
devtools::install_github("hrbrmstr/taucharts")
devtools::install_github("markziemann/Mitch")


# bioc install
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ABSSeq",update=TRUE)
BiocManager::install("edgeR",update=TRUE)
BiocManager::install("DESeq2",update=TRUE)
BiocManager::install("fgsea",update=TRUE)
BiocManager::install("topconfects",update=TRUE)
BiocManager::install("biomaRt",update=TRUE)
BiocManager::install("ExperimentHub",update=TRUE)

BiocManager::install(version='devel')
BiocManager::install("muscat",update=FALSE)
BiocManager::install("scater",update=FALSE)
