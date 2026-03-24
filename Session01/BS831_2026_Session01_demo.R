# BS 831
# In-class demo: Session 1

# check you are in the correct working directory
# use setwd() if not (or in RStudio: Session -> Set Working Directory)
# getwd()


# install packages from CRAN and Bioconductor (if not already installed)

install.packages("BiocManager")
install.packages("tidyverse")

BiocManager::install("Biobase")
BiocManager::install("SummarizedExperiment")



# load packages
library(tidyverse)
library(Biobase)
library(SummarizedExperiment)


# load data
hnsc <- readRDS("C:/Naznin_BU/BU PhD/BU Spring 2026/BS831/HNSC_RNASeq_toy_ES.rds")

# show object
hnsc

# object class, dimensions, etc
class(hnsc)
dim(hnsc)


# ExpressionSet object components

# gene expression data
dim(exprs(hnsc))
exprs(hnsc)[1:5, 1:5]

# phenotype data (i.e. samples / participants, stored in columns)
dim(pData(hnsc))
head(pData(hnsc))
dim(pData(hnsc))

# feature data (i.e. genes, stored in rows)
dim(fData(hnsc))
head(fData(hnsc))


# manipulate ExpressionSet object

# drop second column
pdata1 <- pData(hnsc)[, -2]
head(pdata1)

# using column name
pdata2 <- pData(hnsc)[, !names(pData(hnsc)) %in% "bcr_sample_barcode"]
head(pdata2)

# using tidyverse
pdata3 <- pData(hnsc) |> select(-"bcr_sample_barcode")
head(pdata3)


# SummarizedExperiment

# SummarizedExperiment has additional features compared to ExpressionSet

# convert to SummarizedExperiment
se <- makeSummarizedExperimentFromExpressionSet(hnsc)

class(se)

# show the assay slots (previously exprs())
names(assays(se))

# show colData (previously pData()) and rowData (previously fData())
colData(se)
rowData(se)


# plots

# expression values
hist(assays(se)$exprs)

# log2-transformed expression values
hist(log2(assays(se)$exprs))

# save plots using png() or pdf()
# alternatively use ggsave() from ggplot2
png("../plots/histogram_log2.png", width = 500, height = 500)
hist(log2(assays(se)$exprs))
dev.off()


# additional exercises

# create a slot for extra data
# note log2(0) = -Inf, hence add + 1
assays(se)$log2 <- log2(assays(se)$exprs + 1)

# show object
se

assayNames(se)

# select subset
se_sub <- se[1:100, 1:20]

dim(se_sub)

# all object components are subsetted correctly in the same manner
dim(assays(se_sub)$exprs)
dim(assays(se_sub)$log2)
dim(colData(se_sub))
dim(rowData(se_sub))

