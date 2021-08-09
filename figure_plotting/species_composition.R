library(ggplot2)
library(ggsci)
library(rjson)
library(Cairo)

parameters <- fromJSON(file = "./parameters.json")

fig_dir <- paste0(parameters$publication_dir, "/figs")

pep_species_ratio_6frame <- read.delim(paste(parameters$data_dir, "/accumulated_data/pep_species_ratio.tsv", sep = ""), header=FALSE, stringsAsFactors=FALSE)

trans_composition <- read.delim(paste(parameters$data_dir, "/accumulated_data/reads_species_ratio.tsv", sep = ""), header=FALSE, stringsAsFactors=FALSE)

colnames(pep_species_ratio_6frame) <- c("species", "ratio")
pep_species_ratio_6frame$type <- "Proteogenomics"
sum_ratio <- sum(pep_species_ratio_6frame[,2])
for(i in 1:length(pep_species_ratio_6frame[,2])){
  pep_species_ratio_6frame[i,2] <- pep_species_ratio_6frame[i,2] / sum_ratio
}

colnames(trans_composition) <- c("species", "ratio")
trans_composition$type <- "Transcriptomics"
sum_ratio <- sum(trans_composition[,2])
for(i in 1:length(trans_composition[,2])){
  trans_composition[i,2] <- trans_composition[i,2] / sum_ratio
}

both_composition <- rbind(trans_composition, pep_species_ratio_6frame)

both_composition$species <- unlist(lapply(both_composition$species, FUN = function(x) gsub('\\s', '\n', x)))

Col <- pal_npg("nrc", alpha=0.7)(10)
Col <- c(Col[1], Col[3])
CairoPDF(paste(fig_dir, "/composition_both.pdf", sep = ""), width = 11)
ggplot(both_composition, aes(x=reorder(species, as.numeric(ratio)), y=ratio, fill = type)) +
  geom_bar(stat="identity", position=position_dodge()) +
  coord_flip() +
  labs(y = "Fraction", x = "Bacteria") +
  scale_fill_manual(values = Col) +
	theme(legend.position="bottom",
	      plot.title = element_text(size = 20),
	      axis.text.y = element_text(face = "italic", size = 20),
	      axis.text.x = element_text(size = 25),
	      axis.title = element_text(size = 30),
	      legend.title = element_text(size = 0),
	      legend.text = element_text(size = 30),
	      legend.spacing.x = unit(1, "lines")
	      )
dev.off()
