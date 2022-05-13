#downloading packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")
BiocManager::install(c("affy","affyPLM","sva","AnnotationDbi","hgu133plus2.db"))
library(affy)
library(affyPLM)
library(sva)
library(AnnotationDbi)
library(hgu133plus2.db)
library(bladderbatch)
library(ggplot2)
library(tidyverse)
library(reshape)

#Read in CEL files
data <- ReadAffy(celfile.path = 'samples')

#Normalize CEL files
norm_data <- rma(data, normalize = TRUE)

#Computing Relative Log Expression(not with normalized data)
fit <- fitPLM(data, normalize = TRUE, background = TRUE)
RLE <- RLE(fit)

#Median RLE
RLE_stats <- RLE(fit, type = "stats")


#Computing Normalized Unscaled Standard Errors median (not with normalized data) 
NUSE_stats <- NUSE(fit, type="stats")

#Medians on histogram
hist(RLE_stats["median",], main = "Median RLE stats", xlab = "")
hist(NUSE_stats["median",], main = "Median NUSE stats", xlab = "")

#Combat
annotation <- read.csv('proj_metadata.csv')

batch_correct <- ComBat(exprs(norm_data), 
                        annotation$normalizationcombatbatch, 
                        model.matrix(~normalizationcombatmod, data = annotation))

#csv
write.csv(batch_correct, 'expression_data.csv')

#PCA(must transpose to scale data, then retranspose)
pca <- t(batch_correct) %>% scale() %>%
  t() %>% prcomp(scale=FALSE, center=FALSE)

#Percent variability explained
var_explained <- pca$sdev^2 / sum(pca$sdev^2) *100

#Plotting
plot_data <- as.data.frame(pca$rotation)
plt <- ggplot(plot_data, aes(x=PC1, PC2)) + geom_point() + xlab(paste0("PC1 ", var_explained[1],"%")) + ylab(paste0("PC2 ", var_explained[2],"%")) 

plt
