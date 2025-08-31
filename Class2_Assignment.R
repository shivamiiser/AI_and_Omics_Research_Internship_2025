# function classify_gene() 
classify_gene <- function(logFC, padj) {
  status <- ifelse(padj < 0.05 & logFC > 1, "Upregulated",
                   ifelse(padj < 0.05 & logFC < -1, "Downregulated", "Not_Significant"))
  return(status)
}

#setup directory
setwd("D:/AI_Omics_Internship_2025/Class2")
getwd()

#Set Up the Directorey 
input_dir <- "Raw_Data"
output_dir <- "Results"

# Results directory if it doesn't already exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

files_to_process <- c("DEGs_Data_1.csv", "DEGs_Data_2.csv")


#Process Files Using a for-Loop 



for (file_name in files_to_process) {
  
  cat("Processing file:", file_name, "\n")
  
  input_file_path <- file.path(input_dir, file_name)
  
  dge_data <- read.csv(input_file_path)
  cat("-> File read successfully.\n")
  
  dge_data$padj[is.na(dge_data$padj)] <- 1
  cat("-> Missing padj values replaced with 1.\n")
  
  dge_data$status <- classify_gene(logFC = dge_data$logFC, padj = dge_data$padj)
  cat("-> Gene status classified.\n")
  
  output_file_name <- paste0("classified_", file_name)
  output_file_path <- file.path(output_dir, output_file_name)
  
  write.csv(dge_data, file = output_file_path, row.names = FALSE)
  cat("-> Processed data saved to:", output_file_path, "\n")
  
  cat("\n--- Summary for", file_name, "---\n")
  print(table(dge_data$status))
  cat("\n\n")
}

# Save the entire R workspace
save.image(file = "Shivam_Kumar_Tiwari_Class_2_Assignment.RData")