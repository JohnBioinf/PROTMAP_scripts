library(ggplot2)
library(ggsci)
library(rjson)
library(Cairo)

parameters <- fromJSON(file = "./parameters.json")

fig_dir <- paste0(parameters$publication_dir, "/figs")

data_path <- paste(parameters$data_dir, "/accumulated_data/start_diffs", sep="")
start_diffs <- read.csv(data_path, header=FALSE, stringsAsFactors=FALSE)
colnames(start_diffs) <- c("diff", "k")

start_diffs[,1] <- start_diffs[,1] * -100
Col <- pal_npg("nrc", alpha=0.7)(5)
for(k in c(1, 10)){
	p <- ggplot(start_diffs[start_diffs[,2] > k,], aes(x=diff)) +
		geom_histogram(aes(y=..density..), alpha=0.5,
			       position="identity", bins = 100) + 
		geom_density(alpha=0.5, fill=Col[1]) +
		geom_segment(aes(x = 0, xend = 0, y = -0.01, yend = 0.8),
			     color="red", linetype="dashed") +
		scale_y_continuous(name="Density", limits=c(-0.01, 0.8)) +
		scale_x_continuous(name="Difference in length %",
				   limits=c(-40, 5),
				   breaks=c(-40, -30, -20, -10, 0, 5)) + 
		theme(legend.position="bottom",
		      plot.title=element_text(size=20),
		      axis.text=element_text(size=20),
		      axis.title=element_text(size=30),)
	num_no_diff <- length(start_diffs[(start_diffs[,1] > -5) &
			      (start_diffs[,2] > k),1])
	num_all <- length(start_diffs[start_diffs[,2] > k, 1])
	print(paste(round((num_no_diff/num_all) * 100), "% of proteins can be",
		    "retrieved with 95% of orginal size", sep=""))
	CairoPDF(paste(fig_dir, "/length_compare/size_retrival_", k, ".pdf", sep=""))
	print(p)
	dev.off()
}
