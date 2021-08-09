#!/usr/bin/Rscript

library(ggplot2)
library(ggsci)
library(rjson)
library(Cairo)

parameters <- fromJSON(file = "./parameters.json")

fig_dir <- paste0(parameters$publication_dir, "/figs")

scores_not_annotated_k10_ecoli_path <- paste(parameters$data_dir, "/accumulated_data/scores_not_annotated_k10_ecoli.csv", sep="")
scores_not_annotated_k10_ecoli <- read.csv(scores_not_annotated_k10_ecoli_path, header=FALSE, comment.char="#", stringsAsFactors=T)
colnames(scores_not_annotated_k10_ecoli) <- c("log_score", "category")
scores_not_annotated_k10_ecoli$log_score <- log10(scores_not_annotated_k10_ecoli$log_score)

CairoPDF(paste(fig_dir, "/scores_not_annotated_k10_ecoli.pdf", sep=""), width = 8)
print(ggplot(scores_not_annotated_k10_ecoli, aes(x=log_score, color=category, fill=category)) +
	scale_x_continuous("log(e-value)", lim = c(-12,0), breaks = c(-12,-9,-6,-3,0)) +
	scale_y_continuous("Density") +
        geom_density(alpha=0.5) +
        theme(legend.position="right", legend.direction = "vertical",
              axis.text = element_text(size = 25),
              axis.title = element_text(size = 35),
              legend.title = element_text(size = 0),
              legend.text = element_text(size = 45),
              legend.spacing.x = unit(1, "lines"), 
        )
)
dev.off()
q(save="no")

scores_early_starts_path <- paste(parameters$data_dir, "/accumulated_data/scores_early_starts.csv", sep="")
scores_early_starts <- read.csv(scores_early_starts_path, header=FALSE, comment.char="#", stringsAsFactors=T)
colnames(scores_early_starts) <- c("log_score", "category")
scores_early_starts$log_score <- log10(scores_early_starts$log_score)

CairoPDF(paste(fig_dir, "/scores_early_starts.pdf", sep=""), width = 8)
print(ggplot(scores_early_starts, aes(x=log_score, color=category, fill=category)) +
	scale_x_continuous("log(e-value)") +
	scale_y_continuous("Density") +
        geom_density(alpha=0.5) +
        theme(legend.position="bottom", legend.direction = "vertical",
              axis.text = element_text(size = 25),
              axis.title = element_text(size = 35),
              legend.title = element_text(size = 0),
              legend.text = element_text(size = 45),
              legend.spacing.x = unit(1, "lines"), 
        )
)
dev.off()

