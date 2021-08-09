library(ggsci)
library(eulerr)
library(gridExtra)
library(rjson)
library(Cairo)

names_of_sets = c(" full NCBI anno.", " detected 6frame", " detected NCBI anno.")
parameters <- fromJSON(file = "./parameters.json")

fig_dir <- paste0(parameters$publication_dir, "/figs")

sf1 <- paste(parameters$data_dir, "/accumulated_data/all_CDS", sep = "")
set1 <- array(read.table(sf1)[,1])
k <- 1
Col <- pal_npg("nrc", alpha=0.7)(3)
Col <- c(Col[2], Col[1], Col[3])

for(k in c(1,6,10)){
  sf2 <- paste(parameters$data_dir, "/accumulated_data/k_", k, "_6frame", sep = "")
  sf3 <- paste(parameters$data_dir, "/accumulated_data/k_", k, "_proteom", sep = "")
  
  set2 <- array(read.table(sf2)[,1])
  set3 <- array(read.table(sf3)[,1])
  set_of_sets <- list(set1, set2, set3)
  names(set_of_sets) <- names_of_sets
  
  CairoPDF(paste(fig_dir, "/data_base_compare/venn_", k, ".pdf", sep = ""), width = 7, height = 7)
  p <- plot(euler(set_of_sets, shape = "ellipse"), fills = Col,
       quantities = list(cex = 3),
       legend = list(nrow = 3, ncol = 4, pch = 15, col = Col
                     , labels = names_of_sets, cex = 3.5, side = "bottom")
  )
  grid.arrange(p,widths = c(1,8,1), heights = c(1,8,1),nrow = 3, ncol = 3, layout_matrix = rbind(c(NA, NA, NA),c(NA, 1, NA),c(NA, NA, NA)))
  dev.off()
}
