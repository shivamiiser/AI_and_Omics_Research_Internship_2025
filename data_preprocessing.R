# Bioconductor provides R packages for analyzing omics data (genomics, transcriptomics, proteomics etc).

if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")


# Install CRAN packages for data manipulation
install.packages("dplyr")

# Load Required Libraries
library(GEOquery)             # Download GEO datasets (series matrix, raw CEL files)
library(affy)                 # Pre-processing of Affymetrix microarray data (RMA normalization)
library(arrayQualityMetrics)  # QC reports for microarray data
library(dplyr)                # Data manipulation


#### Raw Data ####
untar("D:/bioinformatics internship/GSE21359_RAW.tar", exdir = "Raw_Data/CEL_Files")


# Read CEL files into R as an AffyBatch object
raw_data <- ReadAffy(celfile.path = "Raw_Data/CEL_Files")

raw_data


#### Quality Control (QC) Before Pre-processing ####

arrayQualityMetrics(expressionset = raw_data,
                    outdir = "Results/QC_Raw_Data",
                    force = TRUE,
                    do.logtransform = TRUE)



#### RMA (Robust Multi-array Average) Normalization ####

normalized_data <- rma(raw_data)

# QC after data normalization 
arrayQualityMetrics(expressionset = normalized_data,
                    outdir = "Results/QC_Normalized_Data",
                    force = TRUE)

# Extract normalized expression values into a data frame
processed_data <- as.data.frame(exprs(normalized_data))

dim(processed_data)   # Dimensions: number of probes × number of samples



#### Filter Low-Variance Transcripts (“soft” intensity based filtering) ####


#Calculate median intensity per probe across samples
row_median <- rowMedians(as.matrix(processed_data))
row_median



# Visualize distribution of probe median intensities
hist(row_median,
     breaks = 100,
     freq = FALSE,
     main = "Median Intensity Distribution")

# Set a threshold to remove low variance probes (dataset-specific, adjust accordingly)
threshold <- 5 
abline(v = threshold, col = "black", lwd = 2) 

# Select probes above threshold
indx <- row_median > threshold 
filtered_data <- processed_data[indx, ] 

gse_data <- getGEO("GSE21359", GSEMatrix = TRUE)[[1]]

# NOW, you can successfully use pData() because 'gse_data' exists
phenotype_data <- pData(gse_data)


# Rename filtered expression data with sample metadata
colnames(filtered_data) <- rownames(phenotype_data)

# Overwrite processed data with filtered dataset
processed_data <- filtered_data 



#### Phenotype Data Preparation ####

class(phenotype_data$source_name_ch1) 

# Define experimental groups (Non-Smoker vs Smoker)
groups <- factor(phenotype_data$`smoking status:ch1`,
                 levels = c("non-smoker", "smoker"),
                 labels = c("Non-Smoker", "Smoker"))
class(groups)
levels(groups)

dim(filtered_data)
