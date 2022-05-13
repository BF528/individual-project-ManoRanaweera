library(tidyverse)
BiocManager::install('tximportData')

expr_data <- read.csv('expression_data.csv') %>%
  as_tibble()

#each gene, at least 20% of the gene-expression values must be > log2(15) (filtered tibble)

filter_20 <- function(tibble){
  tibble$sig_vals <- rowSums(tibble[-1] > (log2(15))) #values per row > log2(15)...counts of this condition per row recorded in extra column "sig_vals"
  filtered_tib <- subset(tibble, ((sig_vals / length(select_if(tibble, is.numeric))) >= 0.2)) #at least 15 % of values meet threshold per row
  return(filtered_tib)
}

filtered <- filter_20(expr_data)

#Have a variance significantly different from the median variance of all probe sets using a threshold of p<0.01

filtered$var <- apply(filtered[,-c(1,136)],1,var)
filtered$median_var <- median(filtered$var)

qchisq_thresh <- qchisq(0.01, 133)


chisq_filtered <- function(tibble){
  tibble$chi <- 133*(filtered$var/filtered$median_var)
  filtered_tib <- subset(tibble, chi > qchisq(0.01,133))%>%
    return()
}

chi_filter_tib <- chisq_filtered(filtered)

chi_filter_tib$sd <- apply(chi_filter_tib[,-c(1,136,137,138,139)],1,sd)

chi_filter_tib$mean <- apply(chi_filter_tib[,-c(1,136,137,138,139)],1,mean)

coeff_filtered <- function(tibble){
  tibble$coeff <- tibble$sd / tibble$mean
  filtered_tib <- subset(tibble, coeff > 0.186) %>%
    return()
}

#Last filter for covariance
coeff_filter_tib <- coeff_filtered(chi_filter_tib)

#Hierarchical clustering

clusters <- t(coeff_filter_tib[,-c(1,136,137,138,139, 140, 141, 142)]) %>%
  scale()%>%
  dist()%>%
  hclust()

plot(clusters,labels = FALSE, main = "Dendrogram", sub="")

#Cut the dendrogram such that the samples are divided into 2 clusters
div_clusters <- cutree(clusters, 2)

#Number of samples in each cluster
table(div_clusters)

#Heatmap and metadata

metadata <- read.csv('proj_metadata.csv')

metadata$color[metadata$cit.coloncancermolecularsubtype == 'C3'] <- 'blue'
metadata$color[metadata$cit.coloncancermolecularsubtype != 'C3'] <- 'red'

heatmap <- heatmap(as.matrix(coeff_filter_tib[,-c(1,136,137,138,139,140,141,142)]), ColSideColors = metadata$color)

coeff_filter_tib <- coeff_filter_tib[,-c(1,136,137,138,139,140,141,142)]

#welch t-test
trans_expr_tib <- coeff_filter_tib 

trans_expr_tib$cluster <- div_clusters

cluster1 <- trans_expr_tib[,div_clusters==1]
cluster2 <- trans_expr_tib[,div_clusters==2]

result <- apply(as.matrix(trans_expr_tib), 1, function(x) t.test(x=x[div_clusters==1],y=x[div_clusters==2]))

tstat <- sapply(result, function(x) x$statistic)
pvalue <- sapply(result, function(x) x$p.value)
padj <- p.adjust(pvalue, method = 'fdr')

diff_exp_genes <- data.frame(coeff_filter_tib$X, tstat, pvalue,padj)

sig_exp_genes <- filter(diff_exp_genes, padj < 0.05)
