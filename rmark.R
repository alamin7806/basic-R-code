## learning R with bioconductor
## load package
library(readr)
library(dplyr)
library(magrittr)
library(tximport)
library(DESeq2)
library(tidyverse)

## not a data.frame but tibble

sample_table = read_csv("SraRunTable.txt") %>% 
  select(`Sample Name`, source_name, TREATMENT,
         Cell_Line,  Cell_type, Time_point) %>%
  slice(c(3,4,13,15))
sample_table = read_csv("SraRunTable.txt")
sample_table = select(sample_table,`Sample Name`, source_name, TREATMENT,
                      Cell_Line,  Cell_type, Time_point)
sample_table = unique(sample_table)
sample_table = slice(sample_table, seq(1,16, by =4))

sample_file = paste0(pull(sample_table,`Sample Name` ), '/quant.sf')

names(sample_file) = 
  pull(sample_table,`Sample Name` )

sample_file[1]
gene_map = read_csv("gene_map.csv", col_names = c('enstid', 'ensgid'))

count_data = tximport(files = sample_file,
                      type = "salmon",
                      tx2gene =gene_map,
                      ignoreTxVersion = TRUE)

dim(count_data[['counts']])

DESeq_DataSet = DESeqDataSetFromTximport(txi = count_data,
                                         colData= sample_table,
                                         design =~conditions )
sample_table = as.data.frame(sample_table)
colnames(sample_table)[1] = "sample"
conditions = c('mock_nhbe','infected_nhbe','mock_4549','infected_a549')
conditions = factor(conditions)
sample_table$conditions = conditions
counts(DESeq_DataSet)[1:6,1:4]
count_data$counts[1:6, 1:4]
DESeq_DataSet = estimateSizeFactors(DESeq_DataSet)
normalizationFactors(DESeq_DataSet)

counts(DESeq_DataSet, normalized =TRUE)[1:6,1:4]

boxplot(counts(DESeq_DataSet, normalized =TRUE))
vst = varianceStabilizingTransformation(DESeq_DataSet)
boxplot(assay(vst))
plotPCA(vst, intgroup = 'conditions') + theme_bw()

dds1 = DESeq_DataSet[,1:2]  
dds2 = DESeq_DataSet[,3:4]

as.matrix(counts(dds1, normalized = FALSE)/ normalizationFactors(dds1))[1:6,1:2]
vst1 = varianceStabilizingTransformation(dds1)
plotPCA(vst1, intgroup = 'conditions')

vst2 = varianceStabilizingTransformation(dds2)
plotPCA(vst2, intgroup = 'conditions')

d = assay(vst1)
d = t(d)
d = dist(d)
h = hclust(d)
plot(h)

s = assay(vst)
s = t(s)
s = dist(s)
s = hclust(s)
plot(s)
k = kmeans(t(assay(vst)), centers = 2)

sample_file_nhbe = sample_file[1:2]

count_data_nhbe = tximport(files = sample_file_nhbe,
                           type = "salmon",
                           tx2gene =gene_map,
                           ignoreTxVersion = TRUE)
sample_table_nhbe = sample_table[1:2,]
sample_table_nhbe$conditions = factor(c('mock','infected'),
                                      levels = c('mock','infected'))

dds_nhbe = DESeqDataSetFromTximport(txi = count_data_nhbe,
                                    colData= sample_table_nhbe,
                                    design =~conditions )
