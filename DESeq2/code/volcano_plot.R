##### volcana_plot #####
volcano_plot <- function( condition = condition, org="3H11", lfc=1, svalue=0.005,ptest=F,pvalue=0.05){
  library("calibrate")

  lfc <- lfc
  ## Test specific directories
  results_dirs= paste("output/", org,"/DEG_tables_svalues/", sep = "")
  results_dirp= paste("output/", org,"/DEG_tables_pvalues/", sep = "")

  # s-value test
  if(ptest == F){
    results_dir <- results_dirs
    p.value = svalue
  }

  # p-value test
  if(ptest == T){
    results_dir <- results_dirp
    p.value = pvalue
  }

  # what analysis is being performed
  cat("\t Plotting volcano plot for:", condition, "\n")

  cond0 <- strsplit(paste(condition, sep=""), split = "condition_")[[1]][2]
  cond1 <- strsplit(cond0, split = "_vs_")[[1]][1]
  cond2 <- strsplit(cond0, split = "_vs_")[[1]][2]
  file <- paste(results_dir, condition, ".txt", sep="")

  df <- read.table(file, sep="\t", header = T)
  if(ptest == F){
    df <- df[order(df$svalue),]
  }
  if(ptest == T){
    df <- df[order(df$padj),]
  }

  # volcano plot
  resPlot <- df
  #resPlot$Gene <- rownames(resPlot)
  resPlot$gene_id <- resPlot$gene_id

  if(ptest == F){
    #p.value <- svalue
    #sig.genes <- dim(subset(resPlot, svalue < p.value & abs(log2FoldChange) > lfc))[1]
    sig.genes <- sum(resPlot$svalue < p.value & abs(resPlot$log2FoldChange) > lfc, na.rm=TRUE)
    outliers <- sum(resPlot$baseMean > 0 & is.na(resPlot$pvalue))
    low.counts <- sum(!is.na(resPlot$pvalue) & is.na(resPlot$svalue))
    minim <- resPlot$log2FoldChange[order(resPlot$log2FoldChange, decreasing = F)][1] + (resPlot$log2FoldChange[order(resPlot$log2FoldChange, decreasing = F)][1]/100) * 20
    maxim <- resPlot$log2FoldChange[order(resPlot$log2FoldChange, decreasing = T)][1] + (resPlot$log2FoldChange[order(resPlot$log2FoldChange, decreasing = T)][1]/100) * 20

    deg <- dim(df[which(abs(df$log2FoldChange) > lfc & df$svalue < p.value),])[1]
    up <- dim(df[which(df$log2FoldChange > lfc & df$svalue < p.value),])[1]
    down <- dim(df[which(df$log2FoldChange < (-lfc) & df$svalue < p.value),])[1]

    with(resPlot, plot(log2FoldChange, -log10(svalue),
                       pch=20, main= paste(condition),
                       sub=paste( sig.genes, " significant genes (L2FC >", lfc, "& svalue <", p.value, ")"),
                       # " | outliers:", outliers,
                       # "| low counts:", low.counts ),
                       xlim=c(minim,maxim),
                       col="gray", cex.sub=0.8, col.sub="gray", cex.lab=0.8))
    # Add colored points: red if svalue< p.value, orange of log2FC > lfc, green if both)
    with(subset(resPlot, svalue < p.value ), points(log2FoldChange, -log10(svalue), pch=20, col="#2c7bb6"))
    with(subset(resPlot, abs(log2FoldChange) > lfc), points(log2FoldChange, -log10(svalue), pch=20, col="#fdae61"))
    with(subset(resPlot, svalue < p.value & (abs(log2FoldChange) > lfc)), points(log2FoldChange, -log10(svalue), pch=20, col="#d7191c"))
    # Label points with the textxy function from the calibrate plot
    if(sig.genes > 9){
      with(subset(resPlot, svalue < p.value & abs(log2FoldChange) > lfc)[1:10,], textxy(log2FoldChange, -log10(svalue), labs=gene_id, cex=.6, col=rgb(0,0,0, 0.5)))
    }
    if(sig.genes <= 9){
      with(subset(resPlot, svalue < p.value & abs(log2FoldChange) > lfc), textxy(log2FoldChange, -log10(svalue), labs=gene_id, cex=.6, col=rgb(0,0,0, 0.5)))
    }
    abline(v=c(-lfc, lfc), h=c(-log10(p.value),-log10(p.value)), col="gray", lty=2)
    mtext(text = paste("UP: ", up), line = 0, adj = 1, col="gray" )
    mtext(text = paste("DOWN: ", down), line = 0, adj = 0, col="gray" )
  }


  if(ptest == T){
    ## p-test plot
    #p.value <- svalue
    #sig.genes <- dim(subset(resPlot, svalue < p.value & abs(log2FoldChange) > lfc))[1]
    sig.genes <- sum(resPlot$padj < p.value & abs(resPlot$log2FoldChange) > lfc, na.rm=TRUE)
    outliers <- sum(resPlot$baseMean > 0 & is.na(resPlot$pvalue))
    low.counts <- sum(!is.na(resPlot$pvalue) & is.na(resPlot$padj))
    minim <- resPlot$log2FoldChange[order(resPlot$log2FoldChange, decreasing = F)][1] + (resPlot$log2FoldChange[order(resPlot$log2FoldChange, decreasing = F)][1]/100) * 20
    maxim <- resPlot$log2FoldChange[order(resPlot$log2FoldChange, decreasing = T)][1] + (resPlot$log2FoldChange[order(resPlot$log2FoldChange, decreasing = T)][1]/100) * 20

    deg <- dim(df[which(abs(df$log2FoldChange) > lfc & df$padj < p.value),])[1]
    up <- dim(df[which(df$log2FoldChange > lfc & df$padj < p.value),])[1]
    down <- dim(df[which(df$log2FoldChange < (-lfc) & df$padj < p.value),])[1]

    with(resPlot, plot(log2FoldChange, -log10(padj),
                       pch=20, main= paste(condition),
                       sub=paste( sig.genes, " significant genes (L2FC >", lfc, "& padj <", p.value, ")"),
                       # " | outliers:", outliers,
                       # "| low counts:", low.counts ),
                       xlim=c(minim,maxim),
                       col="gray", cex.sub=0.8, col.sub="gray", cex.lab=0.8))
    # Add colored points: red if padj< p.value, orange of log2FC > lfc, green if both)
    with(subset(resPlot, padj < p.value ), points(log2FoldChange, -log10(padj), pch=20, col="#2c7bb6"))
    with(subset(resPlot, abs(log2FoldChange) > lfc), points(log2FoldChange, -log10(padj), pch=20, col="#fdae61"))
    with(subset(resPlot, padj < p.value & (abs(log2FoldChange) > lfc)), points(log2FoldChange, -log10(padj), pch=20, col="#d7191c"))
    # Label points with the textxy function from the calibrate plot
    if(sig.genes > 9){
      with(subset(resPlot, padj < p.value & abs(log2FoldChange) > lfc)[1:10,], textxy(log2FoldChange, -log10(padj), labs=gene_id, cex=.6, col=rgb(0,0,0, 0.5)))
    }
    if(sig.genes <= 9){
      with(subset(resPlot, padj < p.value & abs(log2FoldChange) > lfc), textxy(log2FoldChange, -log10(padj), labs=gene_id, cex=.6, col=rgb(0,0,0, 0.5)))
    }
    abline(v=c(-lfc, lfc), h=c(-log10(p.value),-log10(p.value)), col="gray", lty=2)
    mtext(text = paste("UP: ", up), line = 0, adj = 1, col="gray" )
    mtext(text = paste("DOWN: ", down), line = 0, adj = 0, col="gray" )
  }

}
