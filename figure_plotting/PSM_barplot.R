library(ggsci)
library(ggplot2)
library(rjson)
library(Cairo)

draw_key_polygon3 <- function(data, params, size) {
  lwd <- min(data$size, min(size) / 4)
  
  grid::rectGrob(
    width = grid::unit(0.5, "npc"),
    height = grid::unit(0.6, "npc"),
    gp = grid::gpar(
      col = data$colour,
      fill = alpha(data$fill, data$alpha),
      lty = data$linetype,
      lwd = lwd * .pt,
      linejoin = "mitre"
    ))
}

parameters <- fromJSON(file = "./parameters.json")

PSM <- c(rep("1" , 4) , rep("6" , 4), rep("10" , 4))
sub_set_names <- c("only 6frame", "both", "only NCBI", "NCBI not called")
sub_sets <- rep(sub_set_names , 3)
full_NCBI <- paste(parameters$data_dir, "/accumulated_data/all_CDS", sep = "")
full_NCBI <- array(read.table(full_NCBI)[,1])
values <- numeric()
for(k in c(1,6,10)){
  frame <- paste(parameters$data_dir, "/accumulated_data/k_", k, "_6frame", sep = "")
  NCBI <- paste(parameters$data_dir, "/accumulated_data/k_", k, "_proteom", sep = "")
  frame <- array(read.table(frame)[,1])
  NCBI <- array(read.table(NCBI)[,1])
  
  only_frame <-length(setdiff(frame, NCBI))
  both <- length(intersect(frame, NCBI))
  only_NCBI <- length(setdiff(NCBI, frame))
  NCBI_not_called <- length(setdiff(setdiff(full_NCBI, frame), NCBI))
  values <- c(values, only_frame, both, only_NCBI, NCBI_not_called)
}

Col <- pal_npg("nrc", alpha=0.7)(10)
Col <- c(Col[8], Col[9], Col[3], Col[2])

data <- data.frame(PSM,sub_sets,values)
data$sub_sets <- factor(data$sub_sets, levels=sub_set_names)
data$PSM <- factor(data$PSM, levels=c("1", "6", "10"))
p <- ggplot(data, aes(fill=sub_sets, y=values, x=PSM)) + 
  geom_bar(position="stack", stat="identity",key_glyph = "polygon3") +
  scale_fill_manual(values = Col) + 
  geom_segment(x=0.5, xend=3.5, y=length(full_NCBI), yend=length(full_NCBI),
               linetype="dashed", color = "black", size=0.5) +
  labs(y = "Proteins", fill = "") + 
  guides(fill=guide_legend(nrow=2, byrow=FALSE)) +
  theme(legend.position="right",
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 35),
        legend.spacing = unit(1.0, "cm"),
        legend.key.size = unit(4, "cm"),
        legend.margin = margin(0,0,0,200),
        legend.text = element_text(size = 35),
        legend.background=element_blank(),
        legend.key = element_rect(colour = "transparent", fill = "white"))+
CairoPDF(paste("../figs/data_base_compare/barplot_all.pdf", sep=""), width = 8)
print(p)
dev.off()
