library(ggplot2)
library(ggsci)
library(latex2exp)
library(viridis)
library(rjson)
library(Cairo)

load_data <- function(selection_file, data_dir){
  scatterplot_mean <- read.csv(paste0(data_dir, "/accumulated_data/",
                                      selection_file, "_scatterplot_s-hat.csv"),
                               header=FALSE, stringsAsFactors=FALSE)
  colnames(scatterplot_mean) <- c("e_val", "mean_top3", "start_codon", "num_unique_psms",
                                  "genomic_context_translated", "genomic_context_category",
                                  "genomic_context_strand", "transcription", "blast_category",
                                  "name", "candidate")
  scatterplot_mean$num_unique_psms_category <-
  	unlist(lapply(scatterplot_mean$num_unique_psms, FUN = function(x)
  				  if(x > 6){">6"}
  				  else if(x > 2){">2"}
  				  else{as.character(x)}))
  
  scatterplot_mean$num_unique_psms_category <-
  	factor(scatterplot_mean$num_unique_psms_category, levels = c("1", "2", ">2", ">6"))
  scatterplot_mean$blast_category <-
  	factor(scatterplot_mean$blast_category, levels = c("Novel", "Hypothetical", "Known"))
  				
  scatterplot_mean$transcription_category <-
  	unlist(lapply(scatterplot_mean$transcription, FUN = function(x)
  				  if(x > 0.8){">80%"}
  				  else if(x > 0.5){">50%"}
  				  else if(x > 0.1){">10%"}
  				  else{"<10%"}))
  scatterplot_mean$genomic_context_translated_category <-
  	unlist(lapply(scatterplot_mean$genomic_context_translated,
  				  FUN = function(x)
  				  if(x > 6){"PSMs > 6"}
  				  else if(x > 3){"PSMs > 3"}
  				  else if(x > 0){"PSMs > 1"}
  				  else{"no PSMs"}))
  scatterplot_mean$e_val <- log10(scatterplot_mean$e_val) * -1
  scatterplot_mean$mean_top3 <- scatterplot_mean$mean_top3
  return(scatterplot_mean)
}

parameters <- fromJSON(file = "./parameters.json")
fig_dir <- paste0(parameters$publication_dir, "/figs")
current_theme <- theme(legend.position="bottom", legend.direction = "horizontal",
					   axis.text = element_text(size = 30),
					   axis.title = element_text(size = 35),
					   axis.title.y = element_text(angle=0,vjust=0.5),
					   legend.title = element_text(size = 0),
					   legend.text = element_text(size = 30),
					   legend.spacing.x = unit(1, "lines"),
					   legend.key = element_rect(colour = "transparent",
												 fill = "white"),
					   plot.margin = unit(c(0.5,1,0.3,0), "cm"))
guide <- guides(color = guide_legend(override.aes = list(size = 10, alpha = 1), nrow=2, byrow=FALSE))
selection_file <- "nov_psm6"
scatterplot_mean <- load_data(selection_file, parameters$data_dir)

Col <- pal_npg("nrc")(10)
category <- "blast_category"
scatterPlot <- ggplot(data=scatterplot_mean,
                      aes(x=e_val, y=mean_top3, color=eval(parse(text=category)))) +
	scale_colour_manual(values = Col) +
	labs(color=category) +
	geom_point(size=1, alpha=0.3) +
	scale_y_continuous(TeX("$ \\hat{s}$"), lim = c(1,12), breaks = seq(0,12,2)) +
	scale_x_continuous("-log10(e-value)", lim = c(0,20), breaks = seq(0,20,4)) +
	guide +
	current_theme 
CairoPDF(paste0(pic_dir,selection_file,  "_scatterplot_", category, ".pdf"), height=9)
print(scatterPlot)
dev.off()

Col <- inferno(4, begin=0.2, end=0.7)
category <- "transcription_category"
scatterPlot <- ggplot(data=scatterplot_mean,
                      aes(x=e_val, y=mean_top3, color=eval(parse(text=category)))) +
	scale_colour_manual(values = Col) +
	labs(color=category) +
	geom_point(size=1, alpha=0.3) +
	scale_y_continuous(TeX("$ \\hat{s}$"), lim = c(1,12), breaks = seq(0,12,2)) +
	scale_x_continuous("-log10(e-value)", lim = c(0,20), breaks = seq(0,20,4)) +
	guide +
	current_theme 
CairoPDF(paste0(pic_dir, selection_file, "_scatterplot_", category, ".pdf"), height=9)
print(scatterPlot)
dev.off()

Col <- inferno(4, begin=0.2, end=0.7)
category <- "num_unique_psms_category"
scatterPlot <- ggplot(data=scatterplot_mean,
                      aes(x=e_val, y=mean_top3, color=eval(parse(text=category)))) +
	scale_colour_manual(values = Col) +
	labs(color=category) +
	geom_point(size=1, alpha=0.3) +
	scale_y_continuous(TeX("$ \\hat{s}$"), lim = c(1,12), breaks = seq(0,12,2)) +
	scale_x_continuous("-log10(e-value)", lim = c(0,20), breaks = seq(0,20,4)) +
	guide +
	current_theme

CairoPDF(paste0(pic_dir, selection_file, "_scatterplot_", category, ".pdf"), height=9)
print(scatterPlot)
dev.off()

Col <- pal_npg("nrc")(10)
category <- "genomic_context_strand"
scatterPlot <- ggplot(data=scatterplot_mean,
                      aes(x=e_val, y=mean_top3, color=eval(parse(text=category)))) +
	scale_colour_manual(values = Col, labels = c("no\noverlap", "overlap\nother strand", "overlap\nsame strand")) +
	labs(color=category) +
	geom_point(size=1, alpha=0.3) +
	scale_y_continuous(TeX("$ \\hat{s}$"), lim = c(1,12), breaks = seq(0,12,2)) +
	scale_x_continuous("-log10(e-value)", lim = c(0,20), breaks = seq(0,20,4)) +
	guide +
	current_theme
CairoPDF(paste0(pic_dir, selection_file, "_scatterplot_", category, ".pdf"), height=9)
print(scatterPlot)
dev.off()

Col <- inferno(4, begin=0.2, end=0.7)
category <- "genomic_context_translated_category"
scatterPlot <- ggplot(data=scatterplot_mean,
                      aes(x=e_val, y=mean_top3, color=eval(parse(text=category)))) +
	scale_colour_manual(values = Col) +
	labs(color=category) +
	geom_point(size=1, alpha=0.3) +
	scale_y_continuous(TeX("$ \\hat{s}$"), lim = c(1,12), breaks = seq(0,12,2)) +
	scale_x_continuous("-log10(e-value)", lim = c(0,20), breaks = seq(0,20,4)) +
	guide +
	current_theme
CairoPDF(paste0(pic_dir, selection_file, "_scatterplot_", category, ".pdf"), height=9)
print(scatterPlot)
dev.off()
