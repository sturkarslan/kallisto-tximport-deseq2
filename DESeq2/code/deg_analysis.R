##### make_deg: DESeq analysis for all conditions #####

deg_analysis <- function(dds = dds,
                         condition = condition,
                         lfc = 1,
                         write2file = T,
                         sampleTable = sampleTable,
                         svalue = 0.005,
                         org = c("3H11", "R12"),
                         ptest = FALSE) {
  library("reshape2")
  library("DESeq2")

  cat("Performing DEG Analysis for ", org,  "\n")

  # what analysis is being performed
  cat("\t Running DE analysis for condition:", condition, "\n")
  cond0 <- strsplit(paste(condition, sep = ""), split = "condition_")[[1]][2]
  cond1 <- strsplit(cond0, split = "_vs_")[[1]][1]
  cond2 <- strsplit(cond0, split = "_vs_")[[1]][2]

  cat("\t\t Cond1:", cond1, " Cond2:", cond2, "\n")
  cat("\t\t - Running ", org, "analysis\n\n")
  #dds.dv1 <- dds

  # relevel to set cond2 as our reference sample
  cat("... Releveling on ", cond2, "\n")
  dds$condition <- relevel(dds$condition, ref = cond2)

  cat("... Building DESEq object \n")
  dds <- DESeq(dds)

  cat("... Building vsd object \n")
  vsd <- vst(dds, blind = F)

  # get results for s-value
  #resLFC <- results(dds, name = condition, alpha = 0.05)
  # # get results and use lfc shrinkage for visualization
  resLFC <- lfcShrink(
    dds = dds,
    coef = condition,
    type = "apeglm",
    lfcThreshold = lfc
  )
  summary(resLFC)

  # order results and write into a file
  res.ordered <- resLFC[order(resLFC$svalue), ]
  res.ordered <- as_tibble(res.ordered) %>%
    mutate(gene_id = row.names(res.ordered)) %>%
    relocate(gene_id)

  cat("Now writing s-value results to a file\n")
  outputfile = paste("output/",
                     org,
                     "/DEG_tables_svalues/",
                     condition,
                     ".txt",
                     sep = "")
  cat(outputfile, "\n")

  # Add Gene Column
  #res.ordered <- cbind(Transcript = rownames(res.ordered), as.data.frame(res.ordered))
  #res.ordered.final <- merge(res.ordered, as.data.frame(gtf.filtered), by.x ="Transcript", by.y="transcript_id", all.x=T, all.y=F)

  # get summary and write into a file
  if (ptest == F) {
    if (write2file) {
      write_delim(as.data.frame(res.ordered),
                  file = outputfile,
                  delim = "\t")
    }
  }



  # --------------- p-value--------------
  # get results for p-value
  resp <- results(dds, alpha = 0.05, lfcThreshold = lfc, contrast = c("condition", cond1, cond2))

  # order results for padj
  resp.ordered <- resp[order(resp$padj), ]
  resp.ordered <- as.tibble(resp.ordered) %>%
    mutate(gene_id = row.names(resp.ordered)) %>%
    relocate(gene_id)

  # PlotMA-Plot
  cat("MA-Plot for p-value \n")
  plotMA(resp)

  # Plot counts for the most significant gene
  cat("Counts most significant gene for p-value \n")
  plotCounts(dds, gene=which.min(resp$padj), intgroup="condition")


  cat("Now writing p-value results to a file\n")
  outputfile2 = paste("output/",
                      org,
                      "/DEG_tables_pvalues/",
                      condition,
                      ".txt",
                      sep = "")
  cat(outputfile2, "\n")

  # Add Gene Column/Volumes/omics4tb2/ENIGMA/Baliga_pH
  #resp.ordered <- cbind(Transcript = rownames(resp.ordered), as.data.frame(resp.ordered))
  #resp.ordered.final <- merge(resp.ordered, as.data.frame(gtf.filtered), by.x ="Transcript", by.y="transcript_id", all.x=T, all.y=F)

  # get summary and write into a file
  if (ptest == T) {
    if (write2file) {
      write_delim(as.data.frame(resp.ordered),
                  file = outputfile2,
                  delim = "\t")
    }
  }
  if (ptest == F) {
    return(as.data.frame(res.ordered))
  }
  if (ptest == T) {
    return(as.data.frame(resp.ordered))
  }

}
