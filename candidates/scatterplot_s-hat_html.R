library(ggplot2)
library(plotly)
library(htmlwidgets)
library(rjson)
library(Cairo)

args <- commandArgs(trailingOnly=TRUE)
if(is.na(args[1])){
	cat("No Arguments\n")
	q("no", 1)	   
}

parameters <- fromJSON(file = "./parameters.json")

selection_file <- paste0("nov_psm", args[1])

html_dir <- paste0(parameters$data_dir, "/candidates/html/" , selection_file)
data_dir <- paste0(parameters$data_dir, "/accumulated_data/")
scatterplot_mean <- read.csv(paste0(data_dir,
				    selection_file, "_scatterplot_s-hat.csv"),
				    header=FALSE, stringsAsFactors=FALSE)

categorys <- c("start_codon", "genomic_context_category",
	       "genomic_context_strand", "blast_category",
	       "num_unique_psms_category", "transcription_category",
	       "genomic_context_translated_category")

colnames(scatterplot_mean) <- c("e_val", "mean_top3", "start_codon", "num_unique_psms",
                                "genomic_context_translated", "genomic_context_category",
                                "genomic_context_strand", "transcription", "blast_category",
                                "name", "candidate")

scatterplot_mean$num_unique_psms_category <- unlist(lapply(scatterplot_mean$num_unique_psms, FUN = function(x)
  if(x > 6){">6"}else if(x > 2){">2"}else{paste0("=", as.character(x))}))
scatterplot_mean$transcription_category <- unlist(lapply(scatterplot_mean$transcription, FUN = function(x)
  if(x > 0.8){">80%"}else if(x > 0.5){">50%"}else if(x > 0.1){">10%"}else{"<10%"}))
scatterplot_mean$transcription <- scatterplot_mean$transcription * 100
scatterplot_mean$genomic_context_translated_category <- unlist(lapply(scatterplot_mean$genomic_context_translated,
                                                             FUN = function(x)
                                                               if(x > 6){">6PSMs"}
                                                             else if(x > 3){">3PSMs"}
                                                             else if(x > 0){">1PSMs"}
                                                             else{"=0PSMs"}))

scatterplot_mean$e_val <- log10(scatterplot_mean$e_val) * -1
scatterplot_mean$mean_top3 <- scatterplot_mean$mean_top3

scatterplot_mean$url <- paste0("./cand_list.html#content_",
                               unlist(lapply(scatterplot_mean$name,
                                             FUN = function(x) strsplit(x, "_")[[1]][3])))
# category <- "genomic_context_translated"
for(category in categorys){
  print(category)
  scatterPlot <- ggplot(data=scatterplot_mean,
                        aes(x=e_val, y=mean_top3, color=eval(parse(text=category)),
                          text = paste("Name: ", name, "\n", "ID: ", candidate, "\n"),
                          customdata = url)) +
        labs(color=category) +
        geom_point(alpha = 0.4) +
        scale_y_continuous("-log10(mean(e-value))", breaks=seq(1,20,1)) +
        scale_x_continuous("-log10(e-value)", breaks=seq(1,18,1))
  png(paste0(html_dir, "/", category, ".png"), width = 1000, height = 500, type="cairo")
  print(scatterPlot)
  dev.off()
  scatterPlotly <- ggplotly(scatterPlot)
  scatterPlotly <- onRender(scatterPlotly, "
           function(el, x) {
           el.on('plotly_click', function(d) {
           var url = d.points[0].customdata;
           window.open(url);
           });
           }
           ")
  saveWidget(ggplotly(scatterPlotly), file = paste(html_dir, "/", category, ".html", sep=""))
}
# 
category <- "blast_category"
scatterplot_prot <- scatterplot_mean[!duplicated(scatterplot_mean$candidate),]

scatterPlot <- ggplot(data=scatterplot_prot,
                      aes(x=transcription, y=mean_top3, color=eval(parse(text=category)),
                          text = paste("Name: ", name, "\n", "ID: ", candidate, "\n"),
                          customdata = url)) +
      labs(color=category) +
      geom_point(alpha = 0.4) +
        scale_y_continuous("-log10(mean(e-value))", breaks=seq(1,16,1)) +
        scale_x_continuous("% translation", breaks=seq(0,110,10))
png(paste0(html_dir, "/transcription_continuous.png"), width = 1000, height = 500, type="cairo")
print(scatterPlot)
dev.off()
scatterPlotly <- ggplotly(scatterPlot)
scatterPlotly <- onRender(scatterPlotly, "
         function(el, x) {
         el.on('plotly_click', function(d) {
         var url = d.points[0].customdata;
         window.open(url);
         });
         }
         ")
saveWidget(ggplotly(scatterPlotly), file = paste(html_dir, "/transcription_continuous.html", sep=""))

scatterPlot <- ggplot(data=scatterplot_prot,
                      aes(x=genomic_context_translated, y=mean_top3, color=eval(parse(text=category)),
                          text = paste("Name: ", name, "\n", "ID: ", candidate, "\n"),
                          customdata = url)) +
      labs(color=category) +
      geom_point(alpha = 0.4) +
        scale_y_continuous("-log10(mean(e-value))", breaks=seq(1,16,1)) +
        scale_x_continuous("#PSMs", breaks=seq(0,110,10))
png(paste0(html_dir, "/translation_continuous.png"), width = 1000, height = 500, type="cairo")
print(scatterPlot)
dev.off()
scatterPlotly <- ggplotly(scatterPlot)
scatterPlotly <- onRender(scatterPlotly, "
         function(el, x) {
         el.on('plotly_click', function(d) {
         var url = d.points[0].customdata;
         window.open(url);
         });
         }
         ")
saveWidget(ggplotly(scatterPlotly), file = paste(html_dir, "/translation_continuous.html", sep=""))
