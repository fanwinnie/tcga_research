library(bigrquery)
library(tidyverse) #dplyr for manipulating dataframe
                   #ggplot for visualizations

#project name as expressed on BigQuery/GCP

project <- "velvety-ring-155521"

#Avg and SD for miRNA expression on the RIGHT side
r_mirna_exp <- "SELECT mirna_id,
  AVG (reads_per_million_miRNA_mapped) as Avg_Exp_R,
count (*) as num_R,
STDDEV(reads_per_million_miRNA_mapped) as SD_Exp_R
FROM [isb-cgc:TCGA_hg38_data_v0.miRNAseq_Expression]
WHERE case_barcode IN
(SELECT ParticipantBarcode
FROM [eng-name-116319:TCGA_Cancer.SKCM_CSV]
WHERE Asymmetry = 'R')
AND 
sample_barcode IN
(SELECT sample_barcode 
FROM [isb-cgc:TCGA_bioclin_v0.Biospecimen]
WHERE sample_type = '06' OR 
sample_type = '01')
GROUP BY mirna_id"

#execute query and import table
r_mirna_exp <- query_exec(r_mirna_exp, project, max_pages = Inf)

#Avg and SD for miRNA expression on LEFT side
l_mirna_exp <- "SELECT mirna_id,
  AVG (reads_per_million_miRNA_mapped) as Avg_Exp_L,
count (*) as num_L,
STDDEV(reads_per_million_miRNA_mapped) as SD_Exp_L
FROM [isb-cgc:TCGA_hg38_data_v0.miRNAseq_Expression]
WHERE case_barcode IN
(SELECT ParticipantBarcode
FROM [eng-name-116319:TCGA_Cancer.SKCM_CSV]
WHERE Asymmetry = 'L')
AND 
sample_barcode IN
(SELECT sample_barcode 
FROM [isb-cgc:TCGA_bioclin_v0.Biospecimen]
WHERE sample_type = '06' OR 
sample_type = '01')
GROUP BY mirna_id"

#execute query and import table
l_mirna_exp <- query_exec(l_mirna_exp, project, max_pages = Inf)

#merge tables
total_miRNA <- merge(l_mirna_exp, r_mirna_exp, by = "mirna_id", all = TRUE)

#calculate zscore
SE_L = (total_miRNA$SD_Exp_L)/sqrt(total_miRNA$num_L)
SE_R = (total_miRNA$SD_Exp_R)/sqrt(total_miRNA$num_R)
        
total_miRNA$zscore <- ((total_miRNA$Avg_Exp_L - total_miRNA$Avg_Exp_R)/
                       (sqrt((SE_L)^2 + (SE_R)^2)))
hist(total_miRNA$zscore) #check distribution

#filter for significant zscores
miRNA_filtered <- filter(total_miRNA, zscore >= 3.0 | zscore <= -3.0)

#Bonferroni adjustment
  #adjust p-score so that every miRNA is tested independently from the rest
  #this will help correct for the large n
total_miRNA$pvalue2sided <- 2*pnorm(-abs(total_miRNA$zscore))
total_miRNA$bonferroni <- 2*pnorm(-abs(total_miRNA$zscore))*1881
