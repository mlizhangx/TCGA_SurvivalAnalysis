# data on Pancreatic cancer
# library(devtools)
# devtools::install_github(repo = "BioinformaticsFMRP/TCGAbiolinks")

library(TCGAbiolinks)
library(survminer)
library(rms)
library(SummarizedExperiment)
library(stringr)
library(pensim)

setwd("/Users/hai/Documents/workstation/Test_TCGA/TCGA_metadata")
# You can define a list of samples to query and download providing relative TCGA barcodes.
metadata <- read.csv("Prostate cancer_TCGA_metadata.csv")


listSamples <- unique(metadata$aliquot_id)

# Query platform Illumina HiSeq with a list of barcode 
query <- GDCquery(project = "TCGA-PRAD", 
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  experimental.strategy = "RNA-Seq",
                  platform = "Illumina HiSeq",
                  file.type = "results",
                  barcode = listSamples, 
                  legacy = TRUE)
# the prostate should be PRAD, pancreatic is PAAD, and CRC(colorectal) is COAD.

# Download a list of barcodes with platform IlluminaHiSeq_RNASeqV2
# if you already download data, skip this step
GDCdownload(query)

# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
# rsem.genes.results as values
BRCARnaseqSE <- GDCprepare(query)


BRCAMatrix <- assay(BRCARnaseqSE,"raw_count") # or BRCAMatrix <- assay(BRCARnaseqSE,"raw_count")


# data clean on clinical info
dim(metadata) # sample id: 135  meta data feature: 22
count_info_index <- grep("counts", metadata$name)
metadata_count <- metadata[count_info_index, ]
dim(metadata_count) # 45 23

metadata_count$sample_id <- as.character(metadata_count$sample_id)
metadata_count$sample_id <- substr(metadata_count$sample_id,1,nchar(metadata_count$sample_id)-1)


# load survival info
survival_info <- read.csv("prad_tcga_clinical_data.csv") #,sep = "\t"
survival_info <- survival_info[,c("Patient.ID", "Sample.ID", "Overall.Survival..Months.", "Overall.Survival.Status","Gleason.score.overall")]


colnames(survival_info) <- c("pat_id", "sample_id",  "Survival.time",  "status","Gleason_score")
dim(survival_info)
survival_info <- survival_info[!(is.na(survival_info$Survival.time)),]
dim(survival_info)

survival_info$Survival <- ifelse(survival_info$status == "1:DECEASED", 1,0)
survival_info$Gleason_score_binary <- ifelse(survival_info$Gleason_score <= 7, "low","high")

metadata_count_w_survival <- merge(metadata_count,survival_info,by="sample_id")
dim(metadata_count_w_survival)
metadata_count_w_survival$Survival.time.month <- metadata_count_w_survival$Survival.time

# gene matrix: BRCAMatrix
dim(BRCAMatrix) # gene: 19947    samples: 45
gene_count_matrix <- as.data.frame(t(BRCAMatrix))
# match the sample_id with metadata
gene_count_matrix$sample_id <- substr(row.names(gene_count_matrix),1,15)

# merge two data frames by sample_id
dt_merge <- merge(metadata_count_w_survival,gene_count_matrix,by="sample_id")


dt_merge <- dt_merge[dt_merge$sample_type == "Primary Tumor",]
dim(dt_merge)

#### load gene from file # read all gene list ####
read_gene_list <- read.csv("CEMnr_gene_list.csv")
CEMnr_gene_list<- read_gene_list[1:20,]$names


####  select cluster_id and run this section ####
cluster_id <- "CEMnr"
gene_list <- CEMnr_gene_list


length(gene_list) # 1000
length(unique(gene_list)) # 1000
# gene unique gene list
unique_gene_list <- unique(gene_list) 
length(unique_gene_list) # 1000


total_gene_list_in_data <-  str_split_fixed(colnames(gene_count_matrix), "\\|", 2)[,1]

# gene match info:
# count how many gene in gene_list belonge total_gene_list_in_data
length(unique_gene_list) # 1000
length(total_gene_list_in_data) -1  # 19947
sum(unique_gene_list %in% total_gene_list_in_data) # 885

print(paste0("Total # of genes in TCGA: ", nrow(BRCAMatrix)))
print(paste0("Total # of genes in select cluter: ", length(gene_list)))
print(paste0("# of genes match with select cluster: ", sum(unique_gene_list %in% total_gene_list_in_data)))


cph.rt.Survival.covariates.cont <- NULL

for (i in unique_gene_list)
  
{
  sel_gene <- grep(paste0("^",i,"$"), total_gene_list_in_data,value = T)
  # if gene is not in the dataset, jump to next loop
  if(identical(sel_gene, character(0)) == 1) {next}
  print(sel_gene)
  # got the index to get rid of mutiple gene greped from list
  sel_gene_index <- grep(paste0("^",i,"$"), total_gene_list_in_data)
  
  # raw gene count data
  gene_count_raw <- dt_merge[,colnames(gene_count_matrix)[sel_gene_index]]
  # log transfer data
  count_info <- log(gene_count_raw+0.01)
  # standardized data
  count_info <- (gene_count_raw - mean(gene_count_raw)) / sd(gene_count_raw)
  
  cph.rt <- cph(Surv(Survival.time.month,Survival) ~ count_info,  data=dt_merge)
 
  cph.rt.Survival.covariates.cont <- rbind(cph.rt.Survival.covariates.cont,c(cph.rt$n,cph.rt$coefficient,sqrt(cph.rt$var),anova(cph.rt)[1,"P"]))
  
}

rownames(cph.rt.Survival.covariates.cont) <- unique_gene_list[unique_gene_list %in% total_gene_list_in_data]

colnames(cph.rt.Survival.covariates.cont)[3:5] <- c("Est","sd","pvalue")

cph.rt.Survival.covariates.cont <- as.data.frame(cph.rt.Survival.covariates.cont)

cph.rt.Survival.covariates.cont[,"HR"] <- exp(cph.rt.Survival.covariates.cont[,"Est"])

cph.rt.Survival.covariates.cont[,"LHR"] <- exp(cph.rt.Survival.covariates.cont[,"Est"]-1.96*cph.rt.Survival.covariates.cont[,"sd"])

cph.rt.Survival.covariates.cont[,"UHR"] <- exp(cph.rt.Survival.covariates.cont[,"Est"]+1.96*cph.rt.Survival.covariates.cont[,"sd"])


# get the significant gene list from CAMC from survival model using continues value
select_gene_list <- rownames(cph.rt.Survival.covariates.cont)

#### put all significant gene list into one model ####
count_info_table <- data.frame(row.names=1:nrow(dt_merge))

for (i in select_gene_list)
  
{
  sel_gene <- grep(paste0("^",i,"$"), total_gene_list_in_data,value = T)
  # if gene is not in the dataset, jump to next loop
  if(identical(sel_gene, character(0)) == 1) {next}
  print(sel_gene)
  # got the index to get rid of mutiple gene greped from list
  sel_gene_index <- grep(paste0("^",i,"$"), total_gene_list_in_data)
  
  # raw gene count data
  gene_count_raw <- dt_merge[,colnames(gene_count_matrix)[sel_gene_index]]
  # log transfer data
  count_info <- log(gene_count_raw+0.01)
  # standardized data
  count_info <- (gene_count_raw - mean(gene_count_raw)) / sd(gene_count_raw)
  count_info_data <- as.data.frame(count_info)
  colnames(count_info_data) <- paste0(i,"_z")
  count_info_table <- cbind(count_info_table, count_info_data)
}


# get all significant gene in one model
surv_merge_info <- cbind(dt_merge[,c("Survival.time.month","Survival","Gleason_score_binary")], count_info_table)




#### KM plot overall ####
pdf("../Graphs/TCGA_metadata/Prostate KM plot on median Composite survival score CEMnr top 20 gene list with Survival.pdf", height = 8, width = 12)
fit <- npsurv(Surv(Survival.time.month,Survival) ~ 1, data = surv_merge_info)
survplot(fit,n.risk=TRUE,xlab="", lwd=2,main="Overall Survival",conf='bands',time.inc = 12,y.n.risk = -0.2)
fit_summary_table <- as.data.frame(summary(fit)$table)
legend(x=13,y=0.2,legend=paste0("Median OS: NR, 95% CI (", round(fit_summary_table[8,1], 2), " months, ","NR) (n=", nrow(surv_merge_info), ")"),lty=1,bty="n",cex=1) 

fit <- npsurv(Surv(Survival.time.month,Survival) ~ Gleason_score_binary, data = surv_merge_info)
survplot(fit,n.risk=TRUE,xlab="", lwd=2,main="Overall Survival",conf='bands',time.inc = 12,y.n.risk = -0.2,label.curves = FALSE)
fit_summary_table <- as.data.frame(summary(fit)$table)
legend(x=13,y=1.06,legend=paste0("Median Survival on gleason score high: ", round(fit_summary_table[1,7], 2), " months, 95% CI (", round(fit_summary_table[1,8], 2), ", ", round(fit_summary_table[1,9], 2), " ) (n=", nrow(surv_merge_info), ")"),lty=1,bty="n",cex=1) 
legend(x=13,y=1.02,legend=paste0("Median Survival on gleason score low: ", round(fit_summary_table[2,7], 2), " months, 95% CI (", round(fit_summary_table[2,8], 2), ", ", round(fit_summary_table[2,9], 2), " ) (n=", nrow(surv_merge_info), ")"),lty=2,bty="n",cex=1) 

dev.off()
