library(ggsci)
library(ggplot2)
library(latex2exp)
library(rjson)
library(Cairo)
library(gridExtra)
library(grid)


grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}


plot_fdr <- function(fdr_prot,
                     min_k_fdr_prot, min_k_picked_fdr_prot,
                     max_k_fdr_prot, max_k_picked_fdr_prot, title){
	fdr_prot$fdr_prot_log <- unlist(lapply(fdr_prot$fdr_prot,
	                                       FUN = function(x) if(x>0){log10(x)} else {NA}))
	fdr_prot$picked_fdr_prot_log <- unlist(lapply(fdr_prot$picked_fdr_prot,
	                                              FUN = function(x) if(x>0){log10(x)} else {NA}))
	
	lm_fit <- lm(fdr_prot_log ~ k, data = fdr_prot[fdr_prot$k < min_k_fdr_prot &
	                                               fdr_prot$k > max_k_fdr_prot, ])
	inter_classic <- lm_fit$coefficients[1]
	sl_classic <- lm_fit$coefficients[2]
	
	lm_fit <- lm(picked_fdr_prot_log ~ k, data = fdr_prot[fdr_prot$k < min_k_picked_fdr_prot &
	                                                      fdr_prot$k > max_k_picked_fdr_prot, ])
	inter_picked <- lm_fit$coefficients[1]
	sl_picked <- lm_fit$coefficients[2]
	ext_up <- 0.5
	ext_down <- 1
	p <- ggplot(data=fdr_prot, aes(x=k)) +
		# classic
		geom_line(data=fdr_prot[fdr_prot$fdr_prot > 0, ],
		          aes(y = fdr_prot_log, color="a"), size=2) +
		# picked
		geom_line(data=fdr_prot[fdr_prot$picked_fdr_prot > 0, ],
		          aes(y = picked_fdr_prot_log, color="b"), size=2) + 
		# classic
		geom_segment(aes(x = max_k_fdr_prot- ext_up, xend = min_k_fdr_prot + ext_down,
		                 y = sl_classic*(max_k_fdr_prot - ext_up)+inter_classic,
		                 yend = sl_classic*(min_k_fdr_prot + ext_down)+inter_classic,
		                 linetype = "c"), color = "black") +
		# picked
		geom_segment(aes(x = max_k_picked_fdr_prot - ext_up, xend = min_k_picked_fdr_prot + ext_down,
		                 y = sl_picked*(max_k_picked_fdr_prot-ext_up)+inter_picked,
		                 yend = sl_picked*(min_k_picked_fdr_prot+ext_down)+inter_picked,
		                 linetype = "d"), color = "black") +
		scale_x_continuous(TeX("$ \\hat{s}$"), breaks = (0:24)/2, limits = c(0, 4)) +
		scale_y_continuous(bquote("FDR"["prot"]~"%"), limits = c(-2.1,2),
		                   breaks = c(-2, -1, 0, 1, 2),
		                   labels = c("0.01", "0.1", "1.0", "10.0", "100.0")) +
		scale_color_manual(name = "FDR type", values = c(Col[1], Col[2]), labels = c("classical", "picked")) +
		scale_linetype_manual(name = "Regression for FDR types",
		                      values=c("dashed", "dotted"),
		                      labels = c("classical", "picked")) +
		theme(legend.position = "bottom", 
		      axis.text = element_text(size = 20),
		      axis.title = element_text(size = 25),
		      legend.direction = "vertical",
		      legend.title.align = "left",
		      legend.key.width = unit(1.5, "cm"),
		      legend.title = element_text(size = 25),
		      legend.text = element_text(size = 20),
		      plot.title = element_text(size = 30)
		) +
		guides(color = guide_legend(title.position="top", title.hjust = 0.5),
		       linetype = guide_legend(title.position="top", title.hjust = 0.5), 
		       size = guide_legend(nrow = 2, byrow = F)) +
		ggtitle(title)
	return(p)
}

parameters <- fromJSON(file = "./parameters.json")
fig_dir <- paste0(parameters$publication_dir, "/figs")

fdr_prot <- read.csv(paste0(parameters$data_dir, "/accumulated_data/fdr_prot_decoy.csv"),
		     header=FALSE, stringsAsFactors=FALSE)
colnames(fdr_prot) <- c("k", "l" , "fdr_prot", "picked_fdr_prot", "type")
fdr_prot_SIHUMI <- fdr_prot[fdr_prot$l == 6 & fdr_prot$type == "SIHUMI", ]  
fdr_prot_ecoli <- fdr_prot[fdr_prot$l == 6 & fdr_prot$type == "ecoli", ]  

Col <- pal_npg("nrc", alpha = 0.8)(10)

min_k_fdr_prot <- 2.5
min_k_picked_fdr_prot <- 2.5
max_k_fdr_prot <- 1
max_k_picked_fdr_prot <- 1
p_SIHUMI <- plot_fdr(fdr_prot_SIHUMI,
              min_k_fdr_prot, min_k_picked_fdr_prot,
              max_k_fdr_prot, max_k_picked_fdr_prot, "SIHUMIx")

min_k_fdr_prot <- 2
min_k_picked_fdr_prot <- 2
max_k_fdr_prot <- 0.75
max_k_picked_fdr_prot <- 0.75
p_ecoli <- plot_fdr(fdr_prot_ecoli,
              min_k_fdr_prot, min_k_picked_fdr_prot,
              max_k_fdr_prot, max_k_picked_fdr_prot, "E. coli")

CairoPDF(paste(fig_dir, "/FDR_prot_decoy.pdf", sep = ""), width = 7, height = 12)
grid_arrange_shared_legend(p_SIHUMI, p_ecoli)
dev.off()
