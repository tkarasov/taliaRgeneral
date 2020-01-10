# This function takes the results object from deseq and plots with a user defined pvalue cutoff
plot_deseq<-function(res, pval = 0.01, lfc = 1){
  tab = data.frame(logFC = res$log2FoldChange, negLogPval = -log10(res$pvalue))
  head(tab)
  #lfc = 2
  signGenes = (abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval))
  par(mar = c(5, 4, 4, 4))
  p1 <- ggplot(data = tab, aes(x = logFC, y = negLogPval)) +
    geom_point(pch = 16, cex = 2) +
    xlab(expression(log[2]~fold~change)) + 
    ylab(expression(-log[10]~pvalue)) +
    geom_point(data = tab[signGenes, ], aes(x = logFC, y = negLogPval),  pch = 16, cex = 2, col = "red") +
    theme_bw()
  # abline(h = -log10(pval), col = "green3", lty = 2) 
  # abline(v = c(-lfc, lfc), col = "blue", lty = 2) 
  # mtext(paste("pval =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1) 
  # mtext(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc), cex = 0.8, line = 0.5)
  return(p1)
  }

