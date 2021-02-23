## Imports kallisto results into DESEQ2 with tximport
import_kallisto2DESeq <- function(org){
  library('tximport')
  library('rhdf5')
  library("DESeq2")
  library('DT')
  library('tictoc')
  library('apeglm')
  library('viridis')
  library('tidyverse')
  library('Cairo')
  library('stringr')
  setwd("/mnt/omics4tb2-serdar/ENIGMA/Baliga_pH/new_analysis")
  
  cat(" - Running kallisto import \n\n")

  ## Get list of abundance files from kallisto results
  cat("\t - Getting list of abundance files from kallisto \n\n")
  my.dir <- paste("../rnaseq_results_kallisto/",org,sep="")
  sample_id <- dir(file.path(my.dir))
 
  ## list of abundance files.
  files <- file.path(my.dir, sample_id, "abundance.h5")
  names(files) <- sample_id
  cat("\t\t - Collected abundance files for ", length(files), "files in ", org, "\n\n")

  ## Collect sample metadata and prepare sample table
  cat("\t - Collecting meta-data informatioon for samples and creating sample table\n\n")
  ## separate meta data info
  meta1 <- as_tibble(sample_id) %>%
    separate(value, into = c("No","pH1","pH2","Rep","Time"), remove = F, sep = "_") %>%
    mutate(path = paste("../rnaseq_results_kallisto/",org,"/",value,"/abundance.h5",sep = "")) %>%
    mutate(condition = paste(pH1,pH2,Time, sep="_")) %>%
    mutate(treatment = if_else(Time == "Tminus4", "control", "treatment"))
    
    meta2 <- rename(meta1, "sample_name" = value)

  ## rebuild_sample table
  sampleTable <- meta2

  ## Print summary of meta-data info
  print(htmltools::tagList(datatable(sampleTable,rownames = F,escape = F)))
  datatable(sampleTable, caption = "Metadata sample table")

  ## Import kallisto abundances
  cat("\t - Importing count values from collected files for ", org,"\n\n")
  txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE, geneIdCol = "target_id")

  ## Print count matrix header
  print(htmltools::tagList(datatable(txi.kallisto$counts[1:6,1:6],rownames = F,escape = F)))
  datatable(txi.kallisto$counts[1:6,1:6], caption = paste("Kallisto count matrix for ", org))

  ## load data into DESEQ2
  cat("\t - Loading data into DESEQ2 and filtering\n\n")
  dds <- DESeqDataSetFromTximport(txi.kallisto, meta2, ~condition)

  ## Return result objects
  return(list(meta_data=meta1, dds = dds, sample_table= meta2))

}
