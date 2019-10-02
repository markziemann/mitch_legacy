# docker run -d -p 8787:8787 -e USER=rstudio -e PASSWORD=mitch -e ROOT=TRUE rocker/rstudio

RUN apt update && \
  apt install -y \
    zlib1g-dev \
    libxml2-dev \

apt install libxml2-dev




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
                 

# github install
devtools::install_github("markziemann/dee2/getDEE2")
devtools::install_github("hrbrmstr/taucharts")
devtools::install_github("markziemann/Mitch")


# bioc install
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ABSSeq")
BiocManager::install("edgeR")
BiocManager::install("DESeq2")
BiocManager::install("fgsea")
BiocManager::install("topconfects")
BiocManager::install("biomaRt")
BiocManager::install("ExperimentHub")
