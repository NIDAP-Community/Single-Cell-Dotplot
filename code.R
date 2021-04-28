FigureS3d_TableS4 <- function({{{Seurat_Object}}},{{{Cluster_Identities}}}) {
  
  library(tidyverse)
  library(Seurat)
  library(cowplot)
  so <- {{{Seurat_Object}}}$value
  
  #Combine 
  so@meta.data %>% mutate(SCT_snn_res.2.4_mod = case_when(as.numeric(SCT_snn_res.2.8) == 32 ~ 21,
                                                          as.numeric(SCT_snn_res.2.8) == 35 ~ 39,
                                                          as.numeric(SCT_snn_res.2.8) == 41 ~ 36,
                                                          as.numeric(SCT_snn_res.2.8) == 42 ~ 40,
                                                          TRUE ~ as.numeric(SCT_snn_res.2.4))) -> so@meta.data
  Idents(so) <- so@meta.data$SCT_snn_res.2.4_mod
  
  markers <- rev(c("Itgam","Ly6g","Adgre1","Cd68","Siglec1","Csf1r","Apoe","Mafb","Nr4a1",
                   "F13a1","Lgals3","Cxcl9","Cxcl10","Sell","Arg1","Mmp13","Mmp12","Cd38","Mrc1",
                   "Nos2","Cd40","Cd86","Cx3cr1","Cd163","Cd3e","Cd8a","Cd3d","Cd3g","Cd4",
                   "Foxp3","Gzmk","Prf1","Il2ra","Ncr1","Ifng","Klrb1c","Ly6c1","Vcan",
                   "Fn1","Ccr2","S100a9","S100a8","Il1b","G0s2","Csf3r","Cst3","Atox1","Nccrp1",
                   "Ccr9","Siglech","Klk1","Cox6a2","Cd79a","Fcmr"))
  
  dp <- DotPlot(so, assay="SCT", features=markers,dot.scale={{{param1}}},cols = c("lightgrey", "blue"))
  subsamp <- dp$data %>% dplyr::filter(features.plot == "Itgam") %>% arrange(id)
  print(subsamp)
  
  clustertab <- {{{Cluster_Identities}}}
  dp$data$id <- as.character(dp$data$id)
  clusname = clustertab$Latest_Clusnames
  names(clusname) = as.character(clustertab$Clusnum)
  dp$data$name <- clusname[dp$data$id]
  dp$data$name <- factor(dp$data$name,levels=rev(c("NK","CD8","CD4","Tregs","M1_Mac_1","M1_Mac_2","M2_Mac","Mac_1","Mac_2","Mac_3","Mac_4","Mac_5","Mon_1","Mon_2","Mon_3","Mon_4","Mon_5","Mon_6","Mon_7","Mon_8","Mon_9","N_1","N_2","N_3","N_4","N_5","N_7","N_8","DC_1","DC_2","pDC","Unk_1","Unk_2","Unk_3","Unk_4","Unk_5","Unk_6","Unk_7","Unk_8","Unk_9")))
  
  scale.func <- switch(EXPR = "radius", size = scale_size, radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  
  plot <- ggplot(data = dp$data, mapping = aes_string(x = "features.plot", y = "name")) + 
    geom_point(mapping = aes_string(size = "pct.exp", color = "avg.exp.scaled")) + 
    scale.func(range = c(0, 4)) + 
    theme_cowplot() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90))
  print(plot)
  
  # Generate Supplementary Table 4: Contingency Table
  so@meta.data$cellclusters <- clusname[as.character(so@meta.data$SCT_snn_res.2.4_mod)]
  table(so@meta.data$cellclusters,so@meta.data$Likely_CellType)
  
  # Because thresholds for Module Scores were estimated, values may slightly differ from that reported in supplementary table 4.
  
  return(RFoundryObject(so))
}

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in
# the code for any table.

install_bioconductor_package <- function(pkg) {
  install.packages(paste0("https://gypsum.palantircloud.com/assets/dyn/bioconductor-packages/", pkg, ".tar.gz"), repos=NULL)
}
