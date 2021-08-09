library(ggplot2)
library(ggsci)
library(rjson)
library(Cairo)

parameters <- fromJSON(file = "./parameters.json")

scores_only_proteom <- read.table(paste(parameters$data_dir, "/accumulated_data/k_", 1, "_scores_only_proteom.txt", sep=""), quote="\"", comment.char="", stringsAsFactors=FALSE)

Col <- pal_npg("nrc", alpha = 0.6)(10)
Col <- c(Col[9], Col[1], Col[3])

for(k in c(1,6,10)){
	print(k)
	scores_proteom <- read.table(
	  paste(parameters$data_dir, "/accumulated_data/k_", k, "_scores_proteom.txt", sep=""), quote="\"",
	  comment.char="", stringsAsFactors=FALSE)
	scores_6frame <- read.table(
	  paste(parameters$data_dir, "/accumulated_data/k_", k, "_scores_6frame.txt", sep=""), quote="\"",
	  comment.char="", stringsAsFactors=FALSE)

	scores_only_6frame <- read.table(
	  paste(parameters$data_dir, "/accumulated_data/k_", k, "_scores_only_6frame.txt", sep=""), quote="\"",
	  comment.char="", stringsAsFactors=FALSE)
	scores_only_6frame_noCDS <- read.table(
	  paste(parameters$data_dir, "/accumulated_data/k_", k, "_scores_only_6frame_noCDS.txt", sep=""), quote="\"",
	  comment.char="", stringsAsFactors=FALSE)
	scores_only_proteom <- read.table(
	  paste(parameters$data_dir, "/accumulated_data/k_", k, "_scores_only_proteom.txt", sep=""), quote="\"",
	  comment.char="", stringsAsFactors=FALSE)

	scores_both <- c(log10(scores_6frame[,1]), log10(scores_proteom[,1]))
	scores_6frame <- log10(scores_6frame[,1])
	scores_proteom <- log10(scores_proteom[,1])
	scores_only_6frame <- log10(scores_only_6frame[,1])
	scores_only_6frame_noCDS <- log10(scores_only_6frame_noCDS[,1])
	scores_only_proteom <- log10(scores_only_proteom[,1])

# 	scor_fact <- factor(c(rep(1, length(scores_both)),
# 			      rep(2, length(scores_only_6frame)),
# 			      rep(3, length(scores_only_proteom)),
# 			      rep(4, length(scores_only_6frame_noCDS))),
# 			    levels = c(1, 2, 3, 4),
# 			    labels = c("both", "only 6frame", "only proteom", "only 6frame noCDS"))
	scor_fact <- factor(c(rep(1, length(scores_both)),
			      rep(2, length(scores_only_6frame_noCDS)),
			      rep(3, length(scores_only_proteom))),
			    levels = c(1, 2, 3),
			    labels = c("detected both", "detected 6frame", "detected NCBI"))

	scores <- data.frame(data_base = scor_fact,
			     log_score = c(scores_both, scores_only_6frame_noCDS, scores_only_proteom))
	# pdf(paste(parameters$publication_dir, "/figs/data_base_compare/k_", k, "_histogram.pdf", sep=""), width = 8, type="cairo")
	CairoPDF(paste(parameters$publication_dir, "/figs/data_base_compare/k_", k, "_histogram.pdf", sep=""), width = 8)
	p <- ggplot(scores, aes(x=log_score, color=data_base, fill=data_base)) +
	        geom_histogram(alpha=0.5, position="identity") +
	        theme(legend.position="bottom", legend.direction = "vertical",
	              axis.text = element_text(size = 25),
	              axis.title = element_text(size = 35),
	              legend.title = element_text(size = 0),
	              legend.text = element_text(size = 45),
	              legend.spacing.x = unit(1, "lines"), 
	        )
	print(p)
	dev.off()

	# pdf(paste(parameters$publication_dir, "/figs/data_base_compare/k_", k, "_density.pdf", sep=""), width = 8, type="cairo")
	CairoPDF(paste(parameters$publication_dir, "/figs/data_base_compare/k_", k, "_density.pdf", sep=""), width = 8, height = 6)
	print(ggplot(scores, aes(x=log_score, color=data_base, fill=data_base)) +
		scale_x_continuous("log(e-value)", lim = c(-20,0), breaks = c(-20,-15,-10,-5,0)) +
		scale_y_continuous("Density") +
	        geom_density(alpha=0.6) +
		  ggtitle(paste0(k, " PSM")) +
  		scale_fill_manual(values = Col) +
  		scale_color_manual(values = Col) +
	        theme(legend.position="none", legend.direction = "vertical",
	              axis.text = element_text(size = 30),
	              axis.title = element_text(size = 35),
	              plot.title = element_text(size = 50, hjust = 0.5),
	              legend.title = element_text(size = 0),
	              legend.text = element_text(size = 45),
	              legend.spacing.x = unit(1, "lines"), 
	        )
	)
	dev.off()
}
