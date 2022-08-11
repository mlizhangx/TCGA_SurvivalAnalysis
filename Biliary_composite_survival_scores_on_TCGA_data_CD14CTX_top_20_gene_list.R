# use top 20 responder vs NR gene list
# Composite survival scores

library(TCGAbiolinks)
library(survminer)
library(rms)
library(SummarizedExperiment)
library(stringr)
library(pensim)

# setwd("/path/to/your/working/directory/")

# You can define a list of samples to query and download providing relative TCGA barcodes.
metadata <- read.csv("./data/biliarycancer_metadata.csv")
listSamples <- unique(metadata$sample_id)

# Query platform Illumina HiSeq with a list of barcode 
query <- GDCquery(project = "TCGA-CHOL", 
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  experimental.strategy = "RNA-Seq",
                  platform = "Illumina HiSeq",
                  file.type = "results",
                  barcode = listSamples, 
                  legacy = TRUE)

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
dim(metadata_count) 

# load survival info
survival_info <- read.csv("./data/TCGA-CHOL.survival.csv") #,sep = "\t"
colnames(survival_info) <- c("sample_id",    "OS",        "X_PATIENT", "OS.time"  )
metadata_count_w_survival <- merge(metadata_count,survival_info,by="sample_id")
dim(metadata_count_w_survival)
metadata_count_w_survival$OS.time.month <- metadata_count_w_survival$OS.time/30.4

# gene matrix: BRCAMatrix
dim(BRCAMatrix) 
gene_count_matrix <- as.data.frame(t(BRCAMatrix))
# match the sample_id with metadata
gene_count_matrix$sample_id <- substr(row.names(gene_count_matrix),1,16)

# merge two data frames by sample_id
dt_merge <- merge(metadata_count_w_survival,gene_count_matrix,by="sample_id")
dt_merge$status <- ifelse(dt_merge$vital_status == "Dead", 2,1)


dt_merge <- dt_merge[dt_merge$sample_type == "Primary Tumor",]
dim(dt_merge)

### load gene from file # read all gene list ####
read_gene_list <- read.csv("./data/CD14CTX_gene_list.csv", header=TRUE)


colnames(read_gene_list) <- c("row_id","gene","pvalue","logfoldchange","padj")

# run next line if you want sort foldchange by negative values and run pdf() code (line 169) accordingly.
# read_gene_list <- read_gene_list[order(read_gene_list$logfoldchange),]

# run next line if you want sort foldchange by positive values and run pdf() code (line 170) accordingly.
read_gene_list <- read_gene_list[order(-read_gene_list$logfoldchange),]

read_gene_list$names <- read_gene_list$gene

# Top gene list used in model
Top_N = 20
print(Top_N)

Responder_vsNR_T03_top_gene_list <- head(read_gene_list$names, n=Top_N)

####  select cluster_id and run this section ####
cluster_id <- "CD14RvNR"

gene_list <- Responder_vsNR_T03_top_gene_list

length(gene_list) 
length(unique(gene_list)) 
# gene unique gene list
unique_gene_list <- unique(gene_list) 
length(unique_gene_list) 

total_gene_list_in_data <-  str_split_fixed(colnames(gene_count_matrix), "\\|", 2)[,1]

# gene match info:
# count how many gene in gene_list belonge total_gene_list_in_data
length(unique_gene_list) 
length(total_gene_list_in_data) -1 
sum(unique_gene_list %in% total_gene_list_in_data) 

print(paste0("Total # of genes in TCGA: ", nrow(BRCAMatrix)))
print(paste0("Total # of genes in select cluster: ", length(gene_list)))
print(paste0("# of genes match with select cluster: ", sum(unique_gene_list %in% total_gene_list_in_data)))


select_gene_list <- unique_gene_list

#### put all significant gene list into one model ####
# raw count normalized
count_info_table <- data.frame(row.names=1:36)

for (i in select_gene_list)
  
{
  sel_gene <- grep(paste0("^",i,"$"), total_gene_list_in_data,value = T)
  # if gene is not in the dataset, jump to next loop
  if(identical(sel_gene, character(0)) == 1) {next}
  #   print(sel_gene)
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

#### KM plot over all ####
pdf(paste0("./figures/Biliary KM plot on median Composite survival score ", cluster_id, " the most positive fold change in the RvNR ", Top_N," gene list plot.pdf"))

fit <- npsurv(Surv(OS.time.month,OS) ~ 1, data = surv_merge_info)
survplot(fit,n.risk=TRUE,xlab="", lwd=2,main="Overall Survival",conf='bands',time.inc = 3,y.n.risk = -0.2)
fit_summary_table <- as.data.frame(summary(fit)$table)
legend(x=13,y=1,legend=paste0("Median OS: ", round(fit_summary_table[7,1], 2), " months, 95% CI (", round(fit_summary_table[8,1], 2), ", ", round(fit_summary_table[9,1], 2), " ) (n=36)"),lty=1,bty="n",cex=1) 


#### test pensim ####
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
# survfit.training <- survfit(surv.training~preds.training.dichot)
survfit.training <- npsurv(surv.training~preds.training.dichot)
survplot(survfit.training,n.risk=TRUE,xlab="", lwd=2,main="Overall Survival",conf="none",col=c("black", "red"),
         lty =c("dashed","solid"),time.inc = 12,label.curves = FALSE,y.n.risk = -0.2)
summary(coxphfit.training)


(p.training <- signif(summary(coxphfit.training)$logtest[3], 2))  #likelihood ratio test
(hr.training <- signif(summary(coxphfit.training)$conf.int[1], 2))
(hr.lower.training <- summary(coxphfit.training)$conf.int[3])
(hr.upper.training <- summary(coxphfit.training)$conf.int[4])


#### provide median OS for each group rather than HR  
fit_summary_table <- summary(survfit.training)$table

legend(x=19,y=1.06,legend=paste0("Top ", Top_N," gene in model",", Log rank p-value = ", signif(p.training, digits=3)),bty="n",cex=1.0)
legend(x=21,y=1.02,legend=paste0("Median OS for CS > median: ", round(median(fit_summary_table[2,7]), 1), " months"),col="red",lty="solid",bty="n",cex=1.0)
legend(x=21,y=0.98,legend=paste0("Median OS for CS < median: ", round(median(fit_summary_table[1,7]), 1), " months"),col="black",lty="dashed",bty="n",cex= 1.0)

dev.off()

