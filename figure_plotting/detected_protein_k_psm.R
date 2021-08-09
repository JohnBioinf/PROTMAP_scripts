library(ggplot2)
library(ggsci)
library(rjson)
library(Cairo)

parameters <- fromJSON(file = "./parameters.json")

fig_dir <- paste0(parameters$publication_dir, "/figs")

Col <- pal_npg("nrc", alpha = 0.6)(10)
Col <- c(Col[1], Col[3])

discrepancy <- read.delim(paste(parameters$data_dir, "/accumulated_data/detected_protein_k_psm.tsv", sep = ""), header=FALSE, stringsAsFactors=FALSE)
colnames(discrepancy) <- c("k", "num", "type")
discrepancy <- discrepancy[grepl("ecoli", discrepancy[,3], fixed = T),]
CairoPDF(paste0(fig_dir, "/detected_protein_k_psm.pdf"), width = 10)
p <- ggplot(data=discrepancy, aes(x=k, y=num, group=type)) +
	geom_line(aes(color=type), size = 2) +
	scale_x_continuous("Cutoff PSM", 1:20) +
	scale_y_continuous(trans = "log", breaks=c(0,20,50,100,200,400,800)) +
	ylab("Number of Proteins") +
	scale_color_manual(name = "Database", values = c(Col[1], Col[2]),
					   labels = c("only 6frame", "only NCBI")) +
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
