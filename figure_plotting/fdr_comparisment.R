#!/usr/bin/Rscript

library("ggplot2")
library("RColorBrewer")

nth <- function(lst, n){
  sapply(lst, `[`, n)
}

# TODO not adopted old data
colors <- c(brewer.pal(9, "Reds")[4:9], brewer.pal(9,"Blues")[4:9])

fdr_comp <- read.delim("data/ecoli/fdr_comp_new.tsv", header=FALSE, stringsAsFactors = FALSE)
colnames(fdr_comp) <- c("fdr_genome", "fdr_decoy", "experiment")
colnames(fdr_comp) <- c("e-value", "fdr_genome", "fdr_anno", "fdr_decoy", "experiment")

fdr_comp$experiment <- nth(strsplit(fdr_comp$experiment, split = "_"), 3)
fdr_comp$fdr_genome <- fdr_comp$fdr_genome * 1.41

fdr_dup <- fdr_comp
fdr_dup$fdr_anno <- fdr_dup$fdr_genome
fdr_dup$experiment <- paste("fdr_genome", fdr_dup$experiment)
fdr_comp$experiment <- paste("fdr_anno", fdr_comp$experiment)
fdr_new_str <- rbind(fdr_dup[,c(3,4,5)], fdr_comp[, c(3,4,5)])

label_names <- c(expression("1  " ~"FDR"["anno"]),
                 expression("13"~"FDR"["anno"]),
                 expression("17"~"FDR"["anno"]),
                 expression("21"~"FDR"["anno"]),
                 expression("5  " ~"FDR"["anno"]),
                 expression("9  " ~"FDR"["anno"]),
                 expression("1  " ~"FDR"["genome"]),
                 expression("13"~"FDR"["genome"]),
                 expression("17"~"FDR"["genome"]),
                 expression("21"~"FDR"["genome"]),
                 expression("5  " ~"FDR"["genome"]),
                 expression("9  " ~"FDR"["genome"]))

pdf(paste("figs/fdr_comparisment.pdf", sep=""), width = 10)
p <- ggplot(data=fdr_new_str, aes(x=fdr_decoy, y=fdr_anno, group=experiment)) +
	geom_line(aes(color=experiment)) +
	geom_segment(aes(x = 0, xend = 70, y = 0, yend = 70)) +
	scale_colour_manual(values=colors, labels = label_names) +
	labs(color = "Experiments") +
	xlab(bquote("FDR"["decoy"]~"%")) +
	ylab(bquote("FDR"["genome / anno"] ~ "%")) +
	theme(legend.position="right", legend.direction = "vertical",
	      axis.text = element_text(size = 20),
	      axis.title = element_text(size = 25),
	      legend.title = element_text(size = 30),
	      legend.text = element_text(size = 20),
	      legend.spacing.x = unit(1, "lines"),
	      legend.text.align = 0
	)
print(p)
dev.off()
