##### create_DEG_matrix #####
create_DEG_matrix <- function(svalue=0.005,org=c("3H11","R12"), lfc=1, ptest=F, pvalue=0.05,dds){
  cat("Creating DEG Matrix for :", org, "\n")

  # counts results_dir for svalues
  results_dirs= paste("output/", org,"/DEG_tables_svalues/", sep = "")
  # counts results_dir for pvalues
  results_dirp= paste("output/", org,"/DEG_tables_pvalues/", sep = "")


  if(ptest==F){
    testvalue <- svalue
    results_dir <- results_dirs
    sample_ids <- dir(file.path(results_dirs))
    #files <- file.path(results_dir, sample_ids, ".txt")

  }
  if(ptest==T){
    testvalue <- pvalue
    results_dir <- results_dirp
    sample_ids <- dir(file.path(results_dirp))
  }

  # get sample names from all directories for Stylophora
  #sample_id <- unique(dds$meta_data$condition)
   conditions_clean <- vector()
   for(sample in sample_ids){
     cond0 <- strsplit(paste(sample, sep=""), split = "condition_")[[1]][2]
     cond1 <- strsplit(cond0, split = "_vs_")[[1]][1]
     cond2 <- sub(".txt", "", strsplit(cond0, split = "_vs_")[[1]][2])

     conditions_clean <- append(conditions_clean, cond1)
     conditions_clean <- append(conditions_clean, cond2)
   }
   conditions_clean <- unique(conditions_clean)

  DEG.matrix <- matrix(nrow = length(conditions_clean), ncol = length(conditions_clean), dimnames = list(c(conditions_clean), c(conditions_clean)))
  #DEG.matrix <- matrix(nrow = length(sample_ids), ncol = length(sample_ids), dimnames = list(c(sample_ids), c(sample_ids)))

  updown <- data.frame()
  deg_list <- list()
  for(sample in sample_ids){
    cat("\t Processing ", sample, "\n")
    cond0 <- strsplit(paste(sample, sep=""), split = "condition_")[[1]][2]
    cond1 <- strsplit(cond0, split = "_vs_")[[1]][1]
    cond2 <- sub(".txt", "", strsplit(cond0, split = "_vs_")[[1]][2])

    cat(paste(results_dir, sample, ".txt",sep=""))

    df <- read.table(paste(results_dir, sample,sep=""), sep="\t", header = T)

    if(ptest == F){
      deg <- dim(df[which(abs(df$log2FoldChange) > lfc & df$svalue < testvalue),])[1]
      up <- dim(df[which(df$log2FoldChange > lfc & df$svalue < testvalue),])[1]
      up_genes <- rownames(df[which(df$log2FoldChange > lfc & df$svalue < testvalue),][1])
      down <- dim(df[which(df$log2FoldChange < (-lfc) & df$svalue < testvalue),])[1]
      down_genes <- rownames(df[which(df$log2FoldChange < (-lfc) & df$svalue < testvalue),][1])
    }

    if(ptest == T){
      deg <- dim(df[which(abs(df$log2FoldChange) > lfc & df$padj < testvalue),])[1]
      up <- dim(df[which(df$log2FoldChange > lfc & df$padj < testvalue),])[1]
      up_genes <- rownames(df[which(df$log2FoldChange > lfc & df$padj < testvalue),][1])
      down <- dim(df[which(df$log2FoldChange < (-lfc) & df$padj < testvalue),])[1]
      down_genes <- rownames(df[which(df$log2FoldChange < (-lfc) & df$padj < testvalue),][1])
    }


    alter = paste("condition_", cond2, "_vs_", cond1, ".txt", sep="")
    # create names for the lost for listing diferentially expressed genes.
    list_name_up = paste(cond1, "_vs_", cond2, "_UP", sep="")
    list_name_down = paste(cond1, "_vs_", cond2, "_DOWN", sep="")

    if(length(up_genes) > 0){
      deg_list[[list_name_up]] <- up_genes
    }
    if(length(down_genes) > 0){
      deg_list[[list_name_down]] <- down_genes
    }

    updown <- rbind(updown, cbind(comparison=sample, upregulated=up, downregulated=down, cond1=cond1, cond2=cond2, alt=alter))
    DEG.matrix[cond1,cond2] <- deg
  }
  return(list(matrix=DEG.matrix, updown=updown,  deg_genes=deg_list))
}
