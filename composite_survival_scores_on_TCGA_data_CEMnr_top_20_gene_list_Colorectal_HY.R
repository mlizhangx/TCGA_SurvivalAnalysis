# data on Pancreatic cancer
# library(devtools)
# devtools::install_github(repo = "BioinformaticsFMRP/TCGAbiolinks")

library(TCGAbiolinks)
library(survminer)
library(rms)
library(SummarizedExperiment)
library(stringr)
library(pensim)

setwd("/path/to/your/working/directory/")
# You can define a list of samples to query and download providing relative TCGA barcodes.
metadata <- read.csv("Colorectal cancer_TCGA_metadata.csv")


listSamples <- unique(metadata$aliquot_id)

# Query platform Illumina HiSeq with a list of barcode 
query <- GDCquery(project = "TCGA-COAD", 
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
count_info_index <- grep("counts", metadata$name)
metadata_count <- metadata[count_info_index, ]

metadata_count$sample_id <- as.character(metadata_count$sample_id)
metadata_count$sample_id <- substr(metadata_count$sample_id,1,nchar(metadata_count$sample_id)-1)


# load survival info
survival_info <- read.csv("coadread_tcga_clinical_data.csv") #,sep = "\t"
survival_info <- survival_info[,c("Patient.ID", "Sample.ID", "Overall.Survival..Months.", "Overall.Survival.Status")]
colnames(survival_info) <- c("pat_id", "sample_id",  "OS.time",  "status")
survival_info$OS <- ifelse(survival_info$status == "1:DECEASED", 2,1)
table(survival_info$OS )

metadata_count_w_survival <- merge(metadata_count,survival_info,by="sample_id")
dim(metadata_count_w_survival)
metadata_count_w_survival$OS.time.month <- metadata_count_w_survival$OS.time

# gene matrix: BRCAMatrix
gene_count_matrix <- as.data.frame(t(BRCAMatrix))
# match the sample_id with metadata
gene_count_matrix$sample_id <- substr(row.names(gene_count_matrix),1,15)

# merge two data frames by sample_id
dt_merge <- merge(metadata_count_w_survival,gene_count_matrix,by="sample_id")


dt_merge <- dt_merge[dt_merge$sample_type == "Primary Tumor",]

#### load gene from file # read all gene list ####
read_gene_list <- read.csv("CD14_CTX_gene_list.csv")


CD14CTX_gene_list<- read_gene_list[1:20,]$names


####  select cluster_id and run this section ####
cluster_id <- "CD14CTX"
gene_list <- CD14CTX_gene_list


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


cph.rt.OS.covariates.cont <- NULL

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
  
  cph.rt <- cph(Surv(OS.time.month,OS) ~ count_info,  data=dt_merge)

  cph.rt.OS.covariates.cont <- rbind(cph.rt.OS.covariates.cont,c(cph.rt$n,cph.rt$coefficient,sqrt(cph.rt$var),anova(cph.rt)[1,"P"]))
  
}

rownames(cph.rt.OS.covariates.cont) <- unique_gene_list[unique_gene_list %in% total_gene_list_in_data]

colnames(cph.rt.OS.covariates.cont)[3:5] <- c("Est","sd","pvalue")

cph.rt.OS.covariates.cont <- as.data.frame(cph.rt.OS.covariates.cont)

cph.rt.OS.covariates.cont[,"HR"] <- exp(cph.rt.OS.covariates.cont[,"Est"])

cph.rt.OS.covariates.cont[,"LHR"] <- exp(cph.rt.OS.covariates.cont[,"Est"]-1.96*cph.rt.OS.covariates.cont[,"sd"])

cph.rt.OS.covariates.cont[,"UHR"] <- exp(cph.rt.OS.covariates.cont[,"Est"]+1.96*cph.rt.OS.covariates.cont[,"sd"])


# get the significant gene list from CAMC from survival model using continues value
select_gene_list <- rownames(cph.rt.OS.covariates.cont)

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
surv_merge_info <- cbind(dt_merge[,c("OS.time.month","OS")], count_info_table)

noquote(colnames(surv_merge_info))
surv_merge_cox <- cph(Surv(OS.time.month,OS) ~  CXCL5_z+AKR1C1_z+CXCL1_z+SLC7A11_z+CCL7_z+ZNF165_z+MAMLD1_z+SLC16A10_z +  
                      SPINK1_z+BAALC_z+RASAL2_z+TM4SF19_z+CCL2_z+DOCK4_z+ATP8B4_z+NCS1_z+SPRY2_z+MET_z,  data=surv_merge_info)


# get score for each patient
coefficient_list <- surv_merge_cox$coefficient
print(length(coefficient_list))


# adjust the number of coefficient from last output. 
surv_merge_info$score <- surv_merge_info[,3]*coefficient_list[1] + 
  surv_merge_info[,4]*coefficient_list[2] + 
  surv_merge_info[,5]*coefficient_list[3] + 
  surv_merge_info[,6]*coefficient_list[4] + 
  surv_merge_info[,7]*coefficient_list[5] + 
  surv_merge_info[,8]*coefficient_list[6] + 
  surv_merge_info[,9]*coefficient_list[7] +
  surv_merge_info[,10]*coefficient_list[8] +
  surv_merge_info[,11]*coefficient_list[9] +
  surv_merge_info[,12]*coefficient_list[10] +
  surv_merge_info[,13]*coefficient_list[11] +
  surv_merge_info[,14]*coefficient_list[12] +
  surv_merge_info[,15]*coefficient_list[13] +
  surv_merge_info[,16]*coefficient_list[14] +
  surv_merge_info[,17]*coefficient_list[15] +
  surv_merge_info[,18]*coefficient_list[16] +
  surv_merge_info[,19]*coefficient_list[17] +
surv_merge_info[,20]*coefficient_list[18]



surv_merge_cox <- cph(Surv(OS.time.month,OS) ~ score,  data=surv_merge_info)

surv_merge_cox_output <- NULL
surv_merge_cox_output <- data.frame(matrix(ncol = 8, nrow = 1))
colnames(surv_merge_cox_output) <- c("No_Event", "Event",      "Est",     "sd",       "pvalue", "HR","LHR","UHR")

surv_merge_cox_output$No_Event <- surv_merge_cox$n[1]
surv_merge_cox_output$Event <- surv_merge_cox$n[2]
surv_merge_cox_output$Est <- round(surv_merge_cox$coefficient,2)
surv_merge_cox_output$sd <- round(sqrt(surv_merge_cox$var)[1],2)
surv_merge_cox_output$pvalue <- anova(surv_merge_cox)[1,"P"]
surv_merge_cox_output$HR <- round(exp(surv_merge_cox$coefficient),2)
surv_merge_cox_output$LHR <- round(exp(surv_merge_cox$coefficient-1.96*surv_merge_cox_output$sd),2)
surv_merge_cox_output$UHR <- round(exp(surv_merge_cox$coefficient+1.96*surv_merge_cox_output$sd),2)

print(surv_merge_cox_output)

survival.cutoff <- surv_cutpoint(surv_merge_info, time = "OS.time.month", event = "OS", "score",
                                 minprop = 0.1, progressbar = TRUE)

survival.cat <- surv_categorize(survival.cutoff)
colnames(survival.cat) <- c("OS.time.month", "OS", "score_cat")
survival.cat$score <- surv_merge_info$score

print(paste0("cutpoint: ", round(survival.cutoff$cutpoint[1,1],4)))
print(summary(as.numeric(survival.cat$score)))


#### KM plot overall ####
pdf("Colorectal KM plot on median Composite survival score CD14CTX top 20 gene list.pdf")
fit <- npsurv(Surv(OS.time.month,OS) ~ 1, data = surv_merge_info)
survplot(fit,n.risk=TRUE,xlab="", lwd=2,main="Overall Survival",conf='bands',time.inc = 3,y.n.risk = -0.2)
fit_summary_table <- as.data.frame(summary(fit)$table)
legend(x=13,y=1,legend=paste0("Median OS: ", round(fit_summary_table[7,1], 2), " months, 95% CI (", round(fit_summary_table[8,1], 2), ", ", round(fit_summary_table[9,1], 2), " ) (n=", nrow(surv_merge_info), ")"),lty=1,bty="n",cex=1) 



#### KM plot overall ####
pdf("../Graphs/TCGA_metadata/Colorectal KM plot on median Composite survival score CD14CTX top 20 gene list.pdf")
fit <- npsurv(Surv(OS.time.month,OS) ~ 1, data = surv_merge_info)
survplot(fit,n.risk=TRUE,xlab="", lwd=2,main="Overall Survival",conf='bands',time.inc = 12,y.n.risk = -0.2)
fit_summary_table <- as.data.frame(summary(fit)$table)
legend(x=13,y=1,legend=paste0("Median OS: ", round(fit_summary_table[7,1], 2), " months, 95% CI (", round(fit_summary_table[8,1], 2), ", NR) ", "(n=", nrow(surv_merge_info), ")"),lty=1,bty="n",cex=1) 


#### KM plot pensim ####
dim(surv_merge_info) # 253  21
surv_merge_info <- surv_merge_info[!(is.na(surv_merge_info$OS.time.month)),]
dim(surv_merge_info) # 251  21
n_gene <- length(surv_merge_info) - 1

surv.training <- Surv(surv_merge_info$OS.time.month, surv_merge_info$OS)

system.time(output <- opt1D(nsim=1, nprocessors=1, setpen="L2", response=surv.training,
                            penalized=surv_merge_info[,c(3:n_gene)], fold=10, positive=FALSE, standardize=TRUE, 
                            minlambda2=0.2, maxlambda2=100))
# L2" (Ridge) penalty

cc <- output[which.max(output[, "cvl"]), -(1:2)]  #coefficients
sum(abs(cc)>0)  #count non-zero coefficients

preds.training <- as.matrix(surv_merge_info[,c(3:n_gene)]) %*% cc
preds.training.median <- median(preds.training)
preds.training.dichot <- ifelse(preds.training > preds.training.median, "high risk", "low risk")
preds.training.dichot <- factor(preds.training.dichot[, 1], levels=c("low risk", "high risk"))

coxphfit.training <- coxph(surv.training~preds.training.dichot)
survfit.training <- npsurv(surv.training~preds.training.dichot)
summary(coxphfit.training)


(p.training <- signif(summary(coxphfit.training)$logtest[3], 2))  #likelihood ratio test
(hr.training <- signif(summary(coxphfit.training)$conf.int[1], 2))
(hr.lower.training <- summary(coxphfit.training)$conf.int[3])
(hr.upper.training <- summary(coxphfit.training)$conf.int[4])

survplot(survfit.training,n.risk=TRUE,xlab="", lwd=2,main="Overall Survival",conf="none",col=c("black", "red"),
         lty =c("solid","dashed"),time.inc = 12,label.curves = FALSE,y.n.risk = -0.2)

#### provide median OS for each group rather than HR  
fit_summary_table <- summary(survfit.training)$table

survdiff(surv.training~preds.training.dichot)
legend(x=35,y=0.2,legend=paste0("Log rank p-value = 0.002"),bty="n",cex=1.0)
legend(x=21,y=0.15,legend=paste0("Median OS for low CS: NR "),col="black",lty="solid",bty="n",cex= 1.0)
legend(x=21,y=0.1,legend=paste0("Median OS for high CS: ", round(median(fit_summary_table[2,7]), 1), " months"),col="red",lty="dashed",bty="n",cex=1.0)

####

dev.off()







