#!/usr/bin/Rscript

library(ggplot2)
library(RColorBrewer)
library(rjson)
library(Cairo)

# TODO the script that provides the data is not adopted ../data_accumulation/fdr_anno.py
# TODO this script is not adopted

parameters <- fromJSON(file = "./parameters.json")

fig_dir <- paste0(parameters$publication_dir, "/figs")

fdr_prot <- read.csv(paste0(parameters$data_dir, "/fdr_prot.csv"), header=FALSE, stringsAsFactors=FALSE)
colors <- pal_npg("nrc", alpha = 0.8)(2)
colnames(fdr_prot) <- c("k", "fdr_prot", "type")
CairoPDF(paste0(fig_dir, "/fdr_prot.pdf"), width = 10)
p <- ggplot(data=fdr_prot, aes(x=k, y=fdr_prot, group=type)) +
	geom_line(aes(color=type), size = 2) +
	scale_x_continuous("Cutoff PSM", 1:20) +
	scale_y_continuous(bquote("FDR"["prot"]~"%")) +
	scale_colour_manual(values=colors) +
	ylab("") +
	theme(legend.position="bottom", legend.direction = "vertical",
	      axis.text = element_text(size = 20),
	      axis.title = element_text(size = 25),
	      legend.title = element_text(size = 30),
	      legend.text = element_text(size = 20),
	      legend.spacing.x = unit(1, "lines"),
	      legend.text.align = 0
	)
print(p)
dev.off()
